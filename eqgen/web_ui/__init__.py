"""Load the web_ui.html template."""
from pathlib import Path


_HTML_CACHE: str | None = None


def load_html() -> str:
    """Read and return the contents of web_ui.html.

    The result is cached after the first read so repeated calls
    (stdlib server serves each request via a fresh handler) don't
    re-read the file on every page load.
    """
    global _HTML_CACHE
    if _HTML_CACHE is not None:
        return _HTML_CACHE
    path = Path(__file__).parent / "web_ui.html"
    _HTML_CACHE = path.read_text(encoding="utf-8")
    return _HTML_CACHE
