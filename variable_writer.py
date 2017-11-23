from io import StringIO
import jinja2


var_template = jinja2.Template(ur'''
\newcommand\{{varname}}[0]{ {{value}}\xspace}
''')


class VariableCollection(object):
    def __init__(self, name, prefix=None, **vars):
        if prefix is None:
            prefix = name
        self.name = name + '.tex'
        self.prefix = prefix
        self.variables = dict(vars)

    def __getitem__(self, key):
        return self.variables[key]

    def __setitem__(self, key, value):
        self.variables[key] = value

    def define(self, name, value, fp=None):
        varname = "%s%s" % (self.prefix, name)
        text = var_template.render(varname=varname, value=value)
        if fp is None:
            return text
        else:
            fp.write(text)

    def write(self):
        with open(self.name, 'w') as f:
            self.render_definitions(f)

    def render_definitions(self, fp):
        fp.write("% begin variables collection from {self.name}".format(self=self))
        for name, value in sorted(self.variables.items()):
            self.define(name, value, fp)
        fp.write("\n\n")
