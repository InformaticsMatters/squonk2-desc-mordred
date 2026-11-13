"""Put ``src`` on the import path so tests can import the job modules the same
way the container runs them (the scripts live flat in ``/code``)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
