#!/usr/bin/env python3
"""
Extract the complete plain text from paper.tex,
converting all content (abstract, sections, subheadings, paragraphs, captions)
into clean, flowing text without LaTeX formatting, equations, or markup,
ready for rewriting by hand.
"""

import re
import os

def parse_balanced_braces(text, open_idx):
    """Given text with text[open_idx] == '{', return content inside and end index."""
    depth = 0
    for i in range(open_idx, len(text)):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return text[open_idx + 1:i], i + 1
    return text[open_idx + 1:], len(text)

def clean_latex(text):
    # Remove comments
    text = re.sub(r'(?<!\\)%.*?\n', '\n', text)

    # Extract title
    title = ""
    title_match = re.search(r'\\title\{', text)
    if title_match:
        title_raw, _ = parse_balanced_braces(text, title_match.end() - 1)
        title = re.sub(r'\\[a-zA-Z]+', '', title_raw).strip()

    # Get body between \begin{document} and \end{document}
    doc_match = re.search(r'\\begin\{document\}(.*?)\\end\{document\}', text, re.DOTALL)
    if doc_match:
        body = doc_match.group(1)
    else:
        body = text

    # Remove display math \[ ... \] and math environments
    body = re.sub(r'\\\[.*?\\\]', '', body, flags=re.DOTALL)
    body = re.sub(r'\\begin\{equation\*?\}.*?\\end\{equation\*?\}', '', body, flags=re.DOTALL)
    body = re.sub(r'\\begin\{align\*?\}.*?\\end\{align\*?\}', '', body, flags=re.DOTALL)
    body = re.sub(r'\\begin\{gather\*?\}.*?\\end\{gather\*?\}', '', body, flags=re.DOTALL)

    # Replace float environments (figures/tables) with their captions
    def replace_floats(src):
        # find figures/tables
        out = []
        last_idx = 0
        pattern = re.compile(r'\\begin\{(figure|table)\*?\}(.*?)\\end\{\1\*?\}', re.DOTALL)
        for m in pattern.finditer(src):
            out.append(src[last_idx:m.start()])
            env_content = m.group(2)
            # find caption with balanced braces
            cap_m = re.search(r'\\caption\{', env_content)
            if cap_m:
                cap_text, _ = parse_balanced_braces(env_content, cap_m.end() - 1)
                env_type = "Figure" if m.group(1) == "figure" else "Table"
                out.append(f"\n\n[{env_type}: {cap_text.strip()}]\n\n")
            last_idx = m.end()
        out.append(src[last_idx:])
        return "".join(out)

    body = replace_floats(body)

    # Convert inline math $ ... $ to readable text
    def clean_inline_math(m):
        math_content = m.group(1)
        math_replacements = [
            (r'\\pm', '±'),
            (r'\\leq|\\le', '≤'),
            (r'\\geq|\\ge', '≥'),
            (r'\\times', '×'),
            (r'\\approx', '≈'),
            (r'\\sim', '~'),
            (r'\\rightarrow|\\to', '→'),
            (r'\\rho', 'ρ'),
            (r'\\sigma', 'σ'),
            (r'\\alpha', 'α'),
            (r'\\beta', 'β'),
            (r'\\gamma', 'γ'),
            (r'\\mu', 'μ'),
            (r'\\AA', 'Å'),
            (r'\\;', ' '),
            (r'\\,', ' '),
            (r'\\quad', ' '),
            (r'\\math[a-zA-Z]+\{([^}]*)\}', r'\1'),
            (r'\\text[a-zA-Z]*\{([^}]*)\}', r'\1'),
            (r'\\', ''),
            (r'[\{\}]', ''),
        ]
        for p, r_str in math_replacements:
            math_content = re.sub(p, r_str, math_content)
        return math_content.strip()

    body = re.sub(r'\$(.*?)\$', clean_inline_math, body)

    # Clean preamble & document structural markup
    body = re.sub(r'\\maketitle', '', body)
    body = re.sub(r'\\bibliography\{.*?\}', '', body)
    body = re.sub(r'\\bibliographystyle\{.*?\}', '', body)

    # Structure Headings
    def replace_section(m):
        sec_name, _ = parse_balanced_braces(body, m.end() - 1)
        return f"\n\n\n# {sec_name.strip()}\n\n"

    def replace_subsection(m):
        sec_name, _ = parse_balanced_braces(body, m.end() - 1)
        return f"\n\n\n## {sec_name.strip()}\n\n"

    def replace_subsubsection(m):
        sec_name, _ = parse_balanced_braces(body, m.end() - 1)
        return f"\n\n\n### {sec_name.strip()}\n\n"

    # Replace sections with balanced braces
    # Work iteratively
    while True:
        m = re.search(r'\\section\*?\{', body)
        if not m:
            break
        sec_name, end_pos = parse_balanced_braces(body, m.end() - 1)
        body = body[:m.start()] + f"\n\n\n# {sec_name.strip()}\n\n" + body[end_pos:]

    while True:
        m = re.search(r'\\subsection\*?\{', body)
        if not m:
            break
        sec_name, end_pos = parse_balanced_braces(body, m.end() - 1)
        body = body[:m.start()] + f"\n\n\n## {sec_name.strip()}\n\n" + body[end_pos:]

    while True:
        m = re.search(r'\\subsubsection\*?\{', body)
        if not m:
            break
        sec_name, end_pos = parse_balanced_braces(body, m.end() - 1)
        body = body[:m.start()] + f"\n\n\n### {sec_name.strip()}\n\n" + body[end_pos:]

    body = re.sub(r'\\paragraph\*?\{(.*?)\}', r'\n\1: ', body)
    body = re.sub(r'\\begin\{abstract\}', r'\n\n# Abstract\n\n', body)
    body = re.sub(r'\\end\{abstract\}', r'\n\n', body)

    # Citations and cross-references: clean to readable brackets
    body = re.sub(r'\\cite[pt]?\{([^\}]+)\}', r'[\1]', body)
    body = re.sub(r'\\ref\{([^\}]+)\}', r'[\1]', body)
    body = re.sub(r'\\label\{[^\}]+\}', '', body)
    body = re.sub(r'\\url\{([^\}]+)\}', r'\1', body)

    # Formatting macros with balanced braces: \textbf{...}, \emph{...}, etc.
    macro_pattern = re.compile(r'\\(textbf|textit|emph|texttt|small|large|Large|textsc)\{')
    while True:
        m = macro_pattern.search(body)
        if not m:
            break
        inner, end_pos = parse_balanced_braces(body, m.end() - 1)
        body = body[:m.start()] + inner + body[end_pos:]

    # LaTeX entity cleanups
    replacements = [
        (r'~', ' '),
        (r'\{\\AA\}', 'Å'),
        (r'\\AA', 'Å'),
        (r'\\%', '%'),
        (r'\\&', '&'),
        (r'\\\$', '$'),
        (r'\\#', '#'),
        (r'\\_', '_'),
        (r'---', '—'),
        (r'--', '–'),
        (r'``', '"'),
        (r"''", '"'),
        (r'`', "'"),
        (r"'", "'"),
        (r'\\noindent', ''),
        (r'\\small', ''),
        (r'\\centering', ''),
        (r'\\raggedright', ''),
        (r'\\quad', ' '),
        (r'\\qquad', '  '),
        (r'\\\\(\[\d+pt\])?', '\n'),
    ]

    for pat, repl in replacements:
        body = re.sub(pat, repl, body)

    # Remove any leftover commands and isolated braces
    body = re.sub(r'\\[a-zA-Z]+', '', body)
    body = re.sub(r'[\{\}]', '', body)

    # Format clean paragraphs
    raw_blocks = body.split('\n\n')
    cleaned_blocks = []
    
    for block in raw_blocks:
        # Collapse intra-block linebreaks into single spaces
        lines = [l.strip() for l in block.splitlines() if l.strip()]
        if not lines:
            continue
        joined = ' '.join(lines)
        # Normalize multiple spaces
        joined = re.sub(r'\s+', ' ', joined).strip()
        if joined:
            cleaned_blocks.append(joined)

    full_output = f"TITLE: {title}\n\n" + '\n\n'.join(cleaned_blocks)
    return full_output.strip()

if __name__ == '__main__':
    tex_path = '/Volumes/2TB/research-molecule-editing/manuscript/submission/paper.tex'
    out_path = '/Volumes/2TB/research-molecule-editing/manuscript/submission/paper_plain_text.txt'
    
    with open(tex_path, 'r', encoding='utf-8') as f:
        tex_content = f.read()

    plain_text = clean_latex(tex_content)

    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(plain_text)
        f.write('\n')

    # Also save a copy in root workspace for convenience
    root_out = '/Volumes/2TB/research-molecule-editing/paper_plain_text.txt'
    with open(root_out, 'w', encoding='utf-8') as f:
        f.write(plain_text)
        f.write('\n')

    words = plain_text.split()
    print(f"Extraction successful!")
    print(f"Output written to:")
    print(f"  1. {out_path}")
    print(f"  2. {root_out}")
    print(f"Total Word Count: {len(words)}")
