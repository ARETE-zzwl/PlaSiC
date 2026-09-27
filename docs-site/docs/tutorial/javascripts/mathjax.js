// Arithmatex wraps body math, while MkDocs' TOC keeps only the TeX text.
// Scan both: excluding .md-nav would leave raw \(...\) in both desktop and mobile TOCs.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: "tex2jax_ignore|mathjax_ignore",
    processHtmlClass: "arithmatex|tex2jax_process"
  }
};
