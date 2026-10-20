"""Evaluation bundle: hidden truth kept physically separate from the observable bundle (FA-00).

The training path never receives this object.  It carries the hidden expression, the witness (B)
truth and the evaluability mask, and it is written to its own artefact whose name says so.  A
counterexample may change these numbers; the model inputs must not move.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .observation import ObservationContractError, array_digest, json_digest, text_digest


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return float('nan')
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


@dataclass(frozen=True)
class EvaluationBundle:
    """Hidden readouts for one section; never handed to a trainer."""

    sample_id: str
    unit_id: str
    pairing_id: str
    hidden_rows: np.ndarray
    hidden_values: np.ndarray
    hidden_genes: tuple[str, ...]
    witness_genes: tuple[str, ...]
    witness_values: np.ndarray
    truth_kind: str = 'synthetic_counts'
    evaluable_mask: np.ndarray | None = None
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        rows = np.asarray(self.hidden_rows, dtype=np.int64)
        hidden = np.asarray(self.hidden_values, dtype=np.float64)
        witness = np.asarray(self.witness_values, dtype=np.float64)
        if rows.ndim != 1 or rows.size != hidden.shape[0]:
            raise ObservationContractError('hidden_rows must index the rows of hidden_values')
        if hidden.ndim != 2 or hidden.shape[1] != len(self.hidden_genes):
            raise ObservationContractError('hidden_values must be (n_hidden, n_hidden_genes)')
        if witness.ndim != 2 or witness.shape[1] != len(self.witness_genes):
            raise ObservationContractError('witness_values must be (n_rows, n_witness_genes)')
        if not np.all(np.isfinite(hidden)) or np.any(hidden < 0):
            raise ObservationContractError('hidden values must be finite non-negative counts')
        if self.evaluable_mask is not None:
            mask = np.asarray(self.evaluable_mask, dtype=bool)
            if mask.shape != (hidden.shape[0],):
                raise ObservationContractError('evaluable_mask must be per hidden row')
        object.__setattr__(self, 'hidden_rows', rows)
        object.__setattr__(self, 'hidden_values', hidden)
        object.__setattr__(self, 'witness_values', witness)
        object.__setattr__(self, 'hidden_genes', tuple(map(str, self.hidden_genes)))
        object.__setattr__(self, 'witness_genes', tuple(map(str, self.witness_genes)))

    def fingerprint(self) -> dict:
        return dict(sample_id=self.sample_id, unit_id=self.unit_id, pairing_id=self.pairing_id,
                    hidden_rows=array_digest(self.hidden_rows),
                    hidden_values=array_digest(self.hidden_values),
                    hidden_genes=text_digest(list(self.hidden_genes)),
                    witness_genes=text_digest(list(self.witness_genes)),
                    witness_values=array_digest(self.witness_values),
                    truth_kind=self.truth_kind)

    def write(self, path) -> dict:
        path = Path(path)
        if 'evaluation' not in path.name and 'hidden' not in path.name:
            raise ObservationContractError(
                'evaluation artefacts must be named ``*evaluation*`` or ``*hidden*`` so they '
                'cannot be confused with an observable bundle')
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, hidden_rows=np.asarray(self.hidden_rows),
                            hidden_values=np.asarray(self.hidden_values),
                            hidden_genes=np.asarray(self.hidden_genes, dtype=object),
                            witness_genes=np.asarray(self.witness_genes, dtype=object),
                            witness_values=np.asarray(self.witness_values),
                            evaluable_mask=(np.asarray(self.evaluable_mask)
                                            if self.evaluable_mask is not None else np.array([])),
                            sample_id=np.array(self.sample_id), unit_id=np.array(self.unit_id),
                            pairing_id=np.array(self.pairing_id),
                            truth_kind=np.array(self.truth_kind))
        return dict(path=str(path), bytes=path.stat().st_size, fingerprint=self.fingerprint())

    # ------------------------------------------------------------- evaluation only
    def metrics(self, imputed_at_hidden: np.ndarray, prediction_rows: np.ndarray,
                predictions: np.ndarray) -> dict:
        """Quantities that are *allowed* to change when hidden truth or witness truth changes."""
        hidden = self.hidden_values
        mask = (np.ones(hidden.shape[0], dtype=bool) if self.evaluable_mask is None
                else np.asarray(self.evaluable_mask, dtype=bool))
        pred_index = {int(row): i for i, row in enumerate(np.asarray(prediction_rows))}
        paired_pred, paired_truth = [], []
        for local, global_row in enumerate(self.hidden_rows):
            if not mask[local] or int(global_row) not in pred_index:
                continue
            paired_pred.append(float(np.mean(imputed_at_hidden[local])))
            paired_truth.append(float(np.mean(hidden[local])))
        witness = self.witness_values
        b_mean = float(np.mean(witness)) if witness.size else float('nan')
        b_quant = (np.quantile(witness, [0.05, 0.5, 0.95]).astype(float).tolist()
                   if witness.size else [])
        return dict(
            n_hidden=int(hidden.shape[0]), n_evaluable=int(mask.sum()),
            hidden_mean=float(np.mean(hidden)) if hidden.size else float('nan'),
            hidden_var=float(np.var(hidden)) if hidden.size else float('nan'),
            hidden_q50=float(np.median(hidden)) if hidden.size else float('nan'),
            hidden_nonzero_fraction=float(np.mean(hidden > 0)) if hidden.size else float('nan'),
            impute_vs_hidden_corr=_correlation(paired_pred, paired_truth),
            n_impute_pairs=len(paired_pred),
            b_witness_mean=b_mean, b_witness_quantiles=b_quant,
            metrics_sha256=json_digest(dict(
                hidden_mean=float(np.mean(hidden)) if hidden.size else None,
                hidden_var=float(np.var(hidden)) if hidden.size else None,
                impute_vs_hidden_corr=_correlation(paired_pred, paired_truth),
                b_witness_mean=b_mean, b_witness_quantiles=b_quant))['sha256'],
        )
