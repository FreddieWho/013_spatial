#!/usr/bin/env python3
"""Convert the frozen markdown rules into standalone UTF-8 HTML.

`http.server` does not send a charset for `.md`, so browsers decode the Chinese
as something else and show mojibake. Shipping `.html` with an explicit charset
declaration fixes it and renders headings and tables instead of raw markdown.

Deliberately small: only the constructs the rules actually use.
"""

from __future__ import annotations

import argparse
import html
import re
import sys
from pathlib import Path

INLINE_CODE = re.compile(r"`([^`]+)`")
BOLD = re.compile(r"\*\*([^*]+)\*\*")
ITALIC = re.compile(r"(?<!\*)\*([^*\n]+)\*(?!\*)")
LINK = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")
HEADING = re.compile(r"^(#{1,6})\s+(.*)$")
HR = re.compile(r"^\s*(-{3,}|\*{3,})\s*$")
UL = re.compile(r"^(\s*)[-*+]\s+(.*)$")
OL = re.compile(r"^(\s*)(\d+)\.[ \t]+(.*)$")
QUOTE = re.compile(r"^\s*>\s?(.*)$")
TASK = re.compile(r"^\[([ xX])\]\s*(.*)$")
ROW = re.compile(r"^\s*\|.*\|\s*$")
SEP = re.compile(r"^\s*\|[\s:|-]+\|\s*$")

CSS = """
:root { color-scheme: light dark; }
body { margin: 0 auto; max-width: 46em; padding: 32px 22px 80px;
  font: 16px/1.7 -apple-system, "Segoe UI", "Noto Sans CJK SC", "Microsoft YaHei",
        "PingFang SC", sans-serif; }
h1 { font-size: 1.6em; border-bottom: 2px solid #8884; padding-bottom: .3em; }
h2 { font-size: 1.28em; margin-top: 1.8em; border-bottom: 1px solid #8883; padding-bottom: .2em; }
h3 { font-size: 1.08em; margin-top: 1.4em; }
table { border-collapse: collapse; margin: 1em 0; width: 100%; font-size: .94em; }
th, td { border: 1px solid #8885; padding: 6px 10px; text-align: left;
  vertical-align: top; }
th { background: #8881; }
code { background: #8881; padding: 1px 5px; border-radius: 3px;
  font-family: ui-monospace, Menlo, Consolas, monospace; font-size: .9em; }
pre { background: #8881; padding: 12px 14px; border-radius: 6px; overflow-x: auto; }
pre code { background: none; padding: 0; }
blockquote { margin: 1em 0; padding: .2em 1em; border-left: 3px solid #8886;
  color: #aaa; }
hr { border: 0; border-top: 1px solid #8884; margin: 2em 0; }
li { margin: .18em 0; }
a { color: #6aa9ff; }
.backlink { display: inline-block; margin-bottom: 1.6em; color: #888; text-decoration: none; }
"""


def inline(text: str) -> str:
    out = html.escape(text, quote=False)
    out = INLINE_CODE.sub(r"<code>\1</code>", out)
    out = BOLD.sub(r"<strong>\1</strong>", out)
    out = ITALIC.sub(r"<em>\1</em>", out)
    out = LINK.sub(r'<a href="\2">\1</a>', out)
    return out


def convert(md: str, title: str) -> str:
    lines = md.split("\n")
    body: list[str] = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]

        # fenced code
        if line.lstrip().startswith("```"):
            i += 1
            buf = []
            while i < n and not lines[i].lstrip().startswith("```"):
                buf.append(lines[i]); i += 1
            i += 1
            body.append("<pre><code>" + html.escape("\n".join(buf)) + "</code></pre>")
            continue

        # table: header row followed by a separator row
        if ROW.match(line) and i + 1 < n and SEP.match(lines[i + 1]):
            head = [c.strip() for c in line.strip().strip("|").split("|")]
            i += 2
            body.append("<table><thead><tr>"
                        + "".join(f"<th>{inline(c)}</th>" for c in head)
                        + "</tr></thead><tbody>")
            while i < n and ROW.match(lines[i]):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                cells += [""] * (len(head) - len(cells))
                body.append("<tr>" + "".join(f"<td>{inline(c)}</td>"
                                              for c in cells[:len(head)]) + "</tr>")
                i += 1
            body.append("</tbody></table>")
            continue

        m = HEADING.match(line)
        if m:
            lvl = len(m.group(1))
            body.append(f"<h{lvl}>{inline(m.group(2))}</h{lvl}>")
            i += 1
            continue

        if HR.match(line):
            body.append("<hr>")
            i += 1
            continue

        m = QUOTE.match(line)
        if m:
            buf = []
            while i < n and (m := QUOTE.match(lines[i])):
                buf.append(m.group(1)); i += 1
            body.append("<blockquote>" + inline(" ".join(buf)) + "</blockquote>")
            continue

        m = UL.match(line) or OL.match(line)
        if m:
            ordered = OL.match(line) is not None
            tag = "ol" if ordered else "ul"
            base_indent = len(m.group(1))
            items: list[str] = []
            while i < n:
                mm = UL.match(lines[i]) if not ordered else OL.match(lines[i])
                if not mm or len(mm.group(1)) < base_indent:
                    break
                content = mm.group(2) if ordered else mm.group(2)
                t = TASK.match(content)
                if t:
                    checked = " checked" if t.group(1).lower() == "x" else ""
                    content = ('<input type="checkbox" disabled' + checked + "> "
                               + t.group(2))
                    items.append(f'<li class="task">{content}</li>')
                else:
                    items.append(f"<li>{inline(content)}</li>")
                i += 1
            pad = " style='margin-left:1.2em'" if base_indent else ""
            body.append(f"<{tag}{pad}>" + "".join(items) + f"</{tag}>")
            continue

        if not line.strip():
            i += 1
            continue

        buf = []
        while i < n and lines[i].strip() and not (
                HEADING.match(lines[i]) or ROW.match(lines[i]) or HR.match(lines[i])
                or UL.match(lines[i]) or OL.match(lines[i]) or QUOTE.match(lines[i])
                or lines[i].lstrip().startswith("```")):
            buf.append(lines[i].strip())
            i += 1
        if buf:
            body.append("<p>" + inline(" ".join(buf)) + "</p>")
        else:
            i += 1

    return (
        "<!DOCTYPE html>\n<html lang=\"zh-CN\">\n<head>\n"
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{html.escape(title)}</title>\n<style>{CSS}</style>\n</head>\n<body>\n"
        '<a class="backlink" href="../viewer/index.html">&#8592; 返回标注</a>\n'
        + "\n".join(body) + "\n</body>\n</html>\n"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    src, dst = Path(args.src), Path(args.dst)
    md = src.read_text(encoding="utf-8")
    title = args.title or (HEADING.match(md.split("\n")[0]).group(2)
                           if HEADING.match(md.split("\n")[0]) else src.stem)
    dst.write_text(convert(md, title), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())