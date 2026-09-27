"""Keep published note URLs (and their section fragments) working after moves."""

import html
import json
import posixpath
from pathlib import Path
from urllib.parse import urljoin

import yaml
from mkdocs.structure.files import File


def on_post_build(config):
    mappings = yaml.safe_load(
        (Path(__file__).resolve().parents[1] / "redirects.yml").read_text(encoding="utf-8")
    )
    for old, new in mappings.items():
        old_file, new_file = (
            File(path, config.docs_dir, config.site_dir, config.use_directory_urls)
            for path in (old, new)
        )
        if not Path(new_file.abs_dest_path).is_file():
            raise ValueError(f"Redirect destination does not exist: {new}")
        output = Path(old_file.abs_dest_path)
        if output.exists():
            raise ValueError(f"Redirect would overwrite an existing page: {old}")
        relative = posixpath.relpath(new_file.url, posixpath.dirname(old_file.dest_uri))
        if new_file.url.endswith("/"):
            relative += "/"
        # File.url is already URL-encoded, including spaces.
        url = relative
        canonical = html.escape(urljoin(config.site_url, new_file.url), quote=True)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            '<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="robots" content="noindex">'
            f'<link rel="canonical" href="{canonical}">'
            '<title>Page moved</title>'
            f'<script>location.replace({json.dumps(url)} + location.search + location.hash);</script>'
            f'<noscript><meta http-equiv="refresh" content="0; url={url}"></noscript>'
            f'</head><body><p>This page has moved. <a href="{url}">Continue to the new location</a>.</p></body></html>',
            encoding="utf-8",
        )
