# timer.py
import time

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None
        self._end_time = None
        self._time = 0.0

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")
        self._end_time = time.perf_counter()
        elapsed_time = self._end_time - self._start_time
        self._time = elapsed_time

    def print(self):
        if self._start_time is None and self._end_time is None:
            raise TimerError(f"Timer did not run, cannot print")
        if self._start_time is not None and self._end_time is None:
            raise TimerError(f"Timer is running, cannot print.")

        self._start_time = None
        print(f"Elapsed time: {self._time:0.4f} seconds")
