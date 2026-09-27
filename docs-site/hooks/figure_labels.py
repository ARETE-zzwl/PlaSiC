"""把图注开头的 "Figure 2.4.1" / "Fig. 5.1.1." 标签包成 <span class="fig-label">。

网页端由 stylesheets/extra.css 渲染成加粗徽标，PDF 端由 build_pdf.py 的
PRINT_CSS 渲染成打印配色，两边共用同一套标记，保证排版一致。
"""

import re

# 只匹配段落开头的图注标签，正文中段引用的 "Figure 2.4.1" 不受影响。
# 允许 <p> 后紧跟 <em>（图注整体用 *...* 斜体书写的页面）。
FIGURE_LABEL_RE = re.compile(
    r"(<p\b[^>]*>(?:\s*<em\b[^>]*>)?\s*)"
    r"((?:Fig\.|Figure|图)\s*\d+(?:[.\-]\d+)*)"
)


def on_page_content(html, page, config, files):
    return FIGURE_LABEL_RE.sub(r'\1<span class="fig-label">\2</span>', html)
