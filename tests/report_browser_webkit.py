"""Optional real-browser smoke runner; uses system GTK3/WebKit2, not package deps.

Start a temporary broadwayd display, then set GDK_BACKEND=broadway and
BROADWAY_DISPLAY before invoking this with the system Python that provides gi:
    python3 tests/report_browser_webkit.py REPORT.html CHECKS.js --screenshot out.png
Only local file URLs are loaded. No package/JavaScript downloads are needed.
"""
import argparse
from pathlib import Path
import sys


def main():
    import gi
    gi.require_version("Gtk", "3.0")
    gi.require_version("WebKit2", "4.0")
    from gi.repository import Gtk, GLib, WebKit2

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", type=Path)
    parser.add_argument("checks", type=Path)
    parser.add_argument("--screenshot", type=Path)
    parser.add_argument("--fragment", default="")
    parser.add_argument("--no-javascript", action="store_true")
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--height", type=int, default=1050)
    args = parser.parse_args()
    if not Gtk.init_check()[0]:
        raise RuntimeError("No GTK display available; start a temporary broadwayd display first")
    # A headless Broadway display has no physical monitor from which to infer DPI.
    Gtk.Settings.get_default().set_property("gtk-xft-dpi", 96 * 1024)
    window = Gtk.Window()
    window.get_screen().set_resolution(96)
    window.set_default_size(args.width, args.height)
    view = WebKit2.WebView.new_with_context(WebKit2.WebContext.new_ephemeral())
    view.get_settings().set_hardware_acceleration_policy(WebKit2.HardwareAccelerationPolicy.NEVER)
    view.get_settings().set_enable_javascript(not args.no_javascript)
    window.add(view)
    window.show_all()
    state = {"exit": 1, "started": False}

    def fail(error):
        print(str(error), file=sys.stderr, flush=True)
        Gtk.main_quit()

    def snapshot_ready(view, result, _):
        try:
            surface = view.get_snapshot_finish(result)
            surface.write_to_png(str(args.screenshot))
            state["exit"] = 0
            Gtk.main_quit()
        except Exception as error:
            fail(error)

    def snapshot():
        view.get_snapshot(WebKit2.SnapshotRegion.VISIBLE, WebKit2.SnapshotOptions.NONE,
                          None, snapshot_ready, None)
        return False

    def checked(view, result, _):
        try:
            value = view.run_javascript_finish(result).get_js_value().to_string()
            print(value, flush=True)
            if args.screenshot:
                # Allow a rendered frame after the interaction scrolls the page.
                GLib.timeout_add(200, snapshot)
            else:
                state["exit"] = 0
                Gtk.main_quit()
        except Exception as error:
            fail(error)

    def loaded(view, event):
        if event == WebKit2.LoadEvent.FINISHED and not state["started"]:
            state["started"] = True
            # Enabling the evaluator after loading does not execute skipped page scripts.
            if args.no_javascript:
                view.get_settings().set_enable_javascript(True)
            view.run_javascript(args.checks.read_text(encoding="utf-8"), None, checked, None)

    def timed_out():
        fail("Browser checks timed out after 30 seconds")
        return False

    view.connect("load-changed", loaded)
    view.load_uri(args.html.resolve().as_uri() + ("#" + args.fragment if args.fragment else ""))
    GLib.timeout_add_seconds(30, timed_out)
    Gtk.main()
    window.destroy()
    return state["exit"]


if __name__ == "__main__":
    sys.exit(main())
