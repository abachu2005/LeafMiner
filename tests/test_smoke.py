"""Smoke tests for Leaf_Cutter — exercises the bits that are safe to test
without running STAR/LeafCutter2/Slurm in CI."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_sample_fixture_present():
    f = ROOT / "tests" / "fixtures" / "sample.SJ.out.tab"
    assert f.exists()
    rows = [l for l in f.read_text().splitlines() if l.strip()]
    assert len(rows) >= 5
    for line in rows[:5]:
        cols = line.split("\t")
        assert len(cols) == 9, f"STAR SJ.out.tab expects 9 columns, got {len(cols)}"


def test_setup_wizard_imports_clean():
    """The CLI wizard should import without side effects on a fresh Python."""
    import ast
    src = (ROOT / "bin" / "leafcutter-setup").read_text()
    ast.parse(src)


def test_webapp_app_imports():
    """Webapp must import without crashing."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("lc2_app", ROOT / "webapp" / "backend" / "main.py")
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert hasattr(module, "app")
    assert module.app.title.lower().startswith("leafcutter")


def test_frontend_html_has_wizard():
    """Sanity-check the bundled UI still includes the wizard markup."""
    html = (ROOT / "webapp" / "frontend" / "index.html").read_text()
    for marker in ("wizardOverlay", "Setup wizard", "Welcome to LeafCutter2"):
        assert marker in html, f"missing marker: {marker}"


def test_no_personal_paths_committed():
    """Guard against hardcoded NetIDs sneaking back in."""
    import re
    pattern = re.compile(r"iis1026|/projects/p52853/iis1026")
    for p in (ROOT / "scripts").rglob("*.py"):
        text = p.read_text(errors="ignore")
        assert not pattern.search(text), f"Personal path in {p}"
    text = (ROOT / "webapp" / "frontend" / "index.html").read_text()
    assert not pattern.search(text), "Personal path in index.html"
