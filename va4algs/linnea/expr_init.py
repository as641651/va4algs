import pathlib
import os
from tempfile import tempdir
this_dir = pathlib.Path(__file__).parent.resolve()

def expr_init(expr_file, expr_dir):
    template_path = os.path.join(this_dir,'template/generate-variants-linnea.py')
    with open(template_path, "r") as file:
        template_string = file.read()

    with open(expr_file, "r") as file:
        expr_string = file.read()

    for s_ in expr_string.split('\n'):
        if 'def' in s_:  
            expr_fn = s_.split('(')[0].split('def')[-1].strip()
            num_params = len(s_.split('(')[-1].split(','))
            break

    d_str = ''
    for i in range(num_params):
        d_str += '{}_'
        d_str += '{}T/'


    inject = {
        'expression_code':expr_string,
        'num_params':num_params,
        'd_str':d_str,
        'code_dir':expr_dir,
        'expression_fn': expr_fn
    }

    out_file = os.path.join(expr_dir,'generate-variants-linnea.py')

    with open(out_file, "wt", encoding='utf-8') as of:
        of.write(template_string.format(**inject))
    
