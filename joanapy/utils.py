import os


def write_readme_file(params, filename):
    # fout = os.path.join(dir_output, filename)
    fo = open(filename, 'w')
    for k, v in params.items():
        fo.write(str(k) + ' = ' + str(v) + '\n\n')
    fo.close()