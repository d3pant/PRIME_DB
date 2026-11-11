def get_release(self):
    # PRIMEDB's download endpoints send neither Last-Modified nor ETag, which the
    # stock HTTP dumper's set_release() needs -- so the release is the dump date.
    # Imports stay inside the function: the loader inlines this function's source
    # into a generated dumper class.
    import datetime

    return datetime.date.today().isoformat()
