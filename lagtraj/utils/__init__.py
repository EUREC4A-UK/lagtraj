import zlib


def dict_to_hash(d):
    this_hash = 0
    for item in sorted(d.items()):
        curr_hash = 1
        for sub_item in item:
            curr_hash = zlib.adler32(bytes(repr(sub_item), "utf-8"), curr_hash)
        this_hash = this_hash ^ curr_hash

    return str(this_hash)


def optional_debugging(with_debugger):
    """
    Optionally catch exceptions and launch ipdb
    """
    if with_debugger:
        import ipdb

        return ipdb.launch_ipdb_on_exception()
    else:

        class NoDebug:
            def __enter__(self):
                pass

            def __exit__(self, *args, **kwargs):
                pass

        return NoDebug()
