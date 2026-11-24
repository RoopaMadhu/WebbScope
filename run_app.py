"""Entry point used by PyInstaller (and `python run_app.py`). `--mcp` starts Claude's MCP connector (webbscope/mcp.py) instead of the app."""
import sys

if __name__ == "__main__":
    if "--mcp" in sys.argv[1:]:                     # before the app's own imports: the connector needs only the standard library, and its stdout is the protocol
        from webbscope.mcp import main
        main([a for a in sys.argv[1:] if a != "--mcp"])
    else:
        from webbscope.launcher import main
        main()
