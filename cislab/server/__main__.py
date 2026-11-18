"""python -m cislab.server [--port 8000]"""
import argparse
import uvicorn

ap = argparse.ArgumentParser()
ap.add_argument("--port", type=int, default=8000)
ap.add_argument("--host", default="127.0.0.1")
a = ap.parse_args()
print(f"\n  Cislunar Lab is starting.  Open  http://{a.host}:{a.port}  in your browser.\n"
      "  (First start compiles the numerical kernels: allow ~20 seconds.)\n")
uvicorn.run("cislab.server.app:app", host=a.host, port=a.port, log_level="warning")
