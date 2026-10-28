"""CPU tool commands; implementations load only when invoked."""

COMMANDS = {
    'preflight': ('preflight', 'Check input structures, protonation and routing without QM'),
    'evaluate': ('evaluate', 'Evaluate predictions against a declared reference'),
    'build-reference': ('build_reference', 'Build the openTECR-derived reference ledger (needs [reference])'),
}


def command_help():
    return 'Tool commands:\n' + '\n'.join(
        '  metag ' + name + '  — ' + description for name, (_, description) in COMMANDS.items())


def run_command(name, argv):
    from importlib import import_module
    module, _ = COMMANDS[name]
    return import_module('metag.tools.' + module).main(argv)
