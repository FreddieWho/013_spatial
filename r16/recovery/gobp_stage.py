"""Full frozen GO BP validation, using per-spot transforms and visible-only prediction."""
from itertools import product

import numpy as np
from numba import njit
from scipy import sparse
from scipy.spatial import cKDTree

from r16.recovery.corrected import NotTestable, fit_from_moments


@njit(cache=True)
def _scores(data, logged, row, ptr, total, detected, groups, offsets, excludes, eoffsets, normalized):
    n = len(total)
    b = len(eoffsets) - 1
    x = np.ones((n, b, 10))
    y = np.zeros((n, b))
    for j in range(b):
        lib = total.copy()
        nd = detected.copy()
        for q in range(eoffsets[j], eoffsets[j + 1]):
            g = excludes[q]
            for k in range(ptr[g], ptr[g + 1]):
                lib[row[k]] -= data[k]
                nd[row[k]] -= 1
        scale = np.ones(n)
        for i in range(n):
            x[i, j, 1] = np.log1p(max(lib[i], 0.))
            x[i, j, 2] = np.log1p(max(nd[i], 0.))
            if normalized:
                scale[i] = 1e4 / max(lib[i], 1.)
        for a in range(8):
            start, stop = offsets[j * 8 + a], offsets[j * 8 + a + 1]
            values = np.zeros(n)
            for q in range(start, stop):
                g = groups[q]
                for k in range(ptr[g], ptr[g + 1]):
                    value = np.log1p(data[k] * scale[row[k]]) if normalized else logged[k]
                    values[row[k]] += value / (stop - start)
            if a == 7:
                y[:, j] = values
            else:
                x[:, j, a + 3] = values
    return x, y


class FeatureEngine:
    def __init__(self, counts, genes, axes, definitions, unavailable=()):
        self.c = sparse.csc_matrix(counts, dtype=float)
        self.c.eliminate_zeros()
        self.c.sort_indices()
        self.log = np.log1p(self.c.data)
        self.total = np.asarray(self.c.sum(1)).ravel()
        self.detected = np.asarray((self.c > 0).sum(1)).ravel().astype(float)
        self.gx = {g: i for i, g in enumerate(genes)}
        self.available = set(genes) - set(unavailable)
        self.axes = list(axes.values())
        if len(self.axes) != 6:
            raise ValueError('Six composition axes required')
        self.definitions = definitions

    def blocks(self, mode, batch_size=96):
        if mode not in {'raw', 'depth_normalized'}:
            raise ValueError(mode)
        ids = list(self.definitions)
        for start in range(0, len(ids), batch_size):
            valid, bad, groups, offsets, exc, eoffsets = [], [], [], [0], [], [0]
            for name in ids[start:start + batch_size]:
                d = self.definitions[name]
                if set(d['input_genes']) & set(d['readout_genes']):
                    raise ValueError('Input/readout overlap')
                excluded = set(d['input_genes'] + d['readout_genes'])
                if not excluded <= self.available:
                    bad.append((name, 'MISSING_PROGRAM_GENES'))
                    continue
                sets = [sorted(set(a) & self.available - excluded) for a in self.axes]
                if any(not a for a in sets):
                    bad.append((name, 'NOT_SEPARABLE'))
                    continue
                valid.append(name)
                for gs in sets + [d['input_genes'], d['readout_genes']]:
                    groups.extend(self.gx[g] for g in gs)
                    offsets.append(len(groups))
                exc.extend(self.gx[g] for g in sorted(excluded))
                eoffsets.append(len(exc))
            if not valid:
                yield [], None, None, bad
                continue
            x, y = _scores(self.c.data, self.log, self.c.indices, self.c.indptr,
                           self.total, self.detected, np.array(groups), np.array(offsets),
                           np.array(exc), np.array(eoffsets), mode == 'depth_normalized')
            yield valid, x, y, bad


def neighbor_operator(train_xy, query_xy, exclude_self=False):
    k = min(8 + int(exclude_self), len(train_xy))
    d, ix = cKDTree(train_xy).query(query_xy, k=k)
    if k == 1:
        d, ix = d[:, None], ix[:, None]
    if exclude_self:
        d, ix = d[:, 1:], ix[:, 1:]
    if not ix.shape[1]:
        raise NotTestable('INSUFFICIENT_VISIBLE_POINTS')
    weights = 1 / np.maximum(d, 1e-6)
    weights /= weights.sum(1, keepdims=True)
    return sparse.csr_matrix((weights.ravel(), (np.repeat(np.arange(len(query_xy)), ix.shape[1]), ix.ravel())),
                             shape=(len(query_xy), len(train_xy)))


def geometry(xy, hidden):
    hidden = np.asarray(hidden, bool)
    vis = ~hidden
    if hidden.sum() == 0 or vis.sum() < 20:
        raise NotTestable('INSUFFICIENT_VISIBLE_OR_HIDDEN')
    cv, ch = xy[vis], xy[hidden]
    return dict(visible=vis, hidden=hidden, nearest=cKDTree(cv).query(ch)[1],
                neighbors_train=neighbor_operator(cv, cv, True),
                neighbors_test=neighbor_operator(cv, ch))


def moments(x, y):
    xb = np.ascontiguousarray(x.transpose(1, 0, 2))
    return (xb.transpose(0, 2, 1) @ xb / len(x),
            (xb.transpose(0, 2, 1) @ y.T[:, :, None])[:, :, 0] / len(x),
            np.mean(y * y, axis=0))


def mask_predictions(x, y, geo):
    # Never form fit moments from full data and subtract hidden moments: that leaks numerically.
    xv, yv = x[geo['visible']], y[geo['visible']]
    xte = xv[geo['nearest']]
    ntr = geo['neighbors_train'] @ xv[:, :, -1]
    nte = geo['neighbors_test'] @ xv[:, :, -1]
    full = np.concatenate([xv, ntr[:, :, None]], axis=2)
    test = np.concatenate([xte, nte[:, :, None]], axis=2)
    gram, rhs, _ = moments(full, yv)
    preds = []
    for k in [9, 10, 11]:
        beta = fit_from_moments(gram[:, :k, :k], rhs[:, :k])
        preds.append(np.einsum('nbk,bk->nb', test[:, :, :k], beta))
    preds += [yv[geo['nearest']], geo['neighbors_test'] @ yv]
    return np.stack(preds, axis=2)


def sign_flip_p(values):
    """Exact two-sided patient test; requires independent, symmetric patient effects."""
    values = np.asarray(values, float)
    if not len(values) or not np.isfinite(values).all():
        return np.nan
    signs = np.array(list(product([-1., 1.], repeat=len(values))))
    observed = abs(values.sum())
    return float(np.mean(np.abs(signs @ values) >= observed - 1e-12))
