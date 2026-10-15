#!/usr/bin/env python3
"""Печатает переполненные строки (Overfull \\hbox) с привязкой к файлу.
Использование: python3 tools/overfull.py build/glava2.log [порог_pt]"""
import re, sys
log = open(sys.argv[1], encoding='utf-8', errors='replace').read()
thr = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
stack = []
i = 0
tok = re.compile(r'\((\./[^\s()]+\.tex)|\)|Overfull \\hbox \(([\d.]+)pt too wide\) (?:in paragraph|detected) at lines? (\d+)')
for m in tok.finditer(log):
    if m.group(1):
        stack.append(m.group(1))
    elif m.group(0) == ')':
        if stack: stack.pop()
    else:
        w = float(m.group(2))
        if w >= thr:
            print(f"{stack[-1] if stack else '?'}:{m.group(3)}  {w:.1f}pt")
