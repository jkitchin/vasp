def monkeypatch_class(cls):
    """Decorator to add methods to a class.

    @monkeypatch_class(Vasp)
    def some_func(self):
        return "you got it."

    """
    def decorator(func):
        setattr(cls, func.__name__, func)
        s = ('\n\nMonkey-patch defined in '
             '{f.func_code.co_filename} '
             'at line {f.func_code.co_firstlineno}')
        if func.__doc__ is None:
            func.__doc__ = ''
            # I am not happy about adding this try here, but some
            # functions seem to be missing the func_code attribute and
            # I am covering this here.
            try:
                func.__doc__ += s.format(f=func)
            except:
                pass
        return func
    return decorator
