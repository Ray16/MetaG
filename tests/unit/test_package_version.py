"""Ensure the source and distribution do not advertise different versions."""
import re
from pathlib import Path
import metag

def test_package_version_matches_build_metadata():
    text=(Path(__file__).resolve().parents[2]/'pyproject.toml').read_text()
    versions=re.findall(r'^version\s*=\s*"([^"]+)"\s*$',text,flags=re.MULTILINE)
    assert versions==[metag.__version__]
