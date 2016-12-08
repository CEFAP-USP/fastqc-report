#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Workflow:
# 1. Check if the folder has been analysed before
#   1.1 Status: checked, converted, reported, compiled, emailed, running, error, completed
# 2. If the sequencer is NextSeq:
#   2.1 Run bcl2fastq to create the FASTQ files
#   2.1.1 Execution:
#       nohup /usr/local/bin/bcl2fastq
#           --runfolder-dir 160225_NB501279_0002_AHTGNYBGXX/
#           --output-dir 160225_NB501279_0002_AHTGNYBGXX_fastq &
# 3. Run FastQC with the files created on output-dir on 2.1
#   3.1 /data/runs/FastQC/FastQC/fastqc --extract -t 8 Undetermined_S0_L00?_R1_001.fastq.gz
# 4. Compile tex with the results on 3.1
#   4.1 pdflatex -output-directory [DIR] tex.tex
# 5. Send email with the PDF created on 4.1
#   5.1 sendmail ...

import argparse
import os
import subprocess
import shutil
import csv
import re
from collections import OrderedDict
from bs4 import BeautifulSoup
import datetime


BCL2FASTQ_PATH = '/usr/local/bin/bcl2fastq'
FASTQC_PATH = '/data/runs/FastQC/FastQC/fastqc'
WORKING_DIR = os.path.dirname(os.path.abspath(__file__))
REPORT_FILE = 'FastQC_report.tex'
REPORTS_PATH = 'FastQC_reports'
STATUS_FILE = 'run_report'
# informações do experimento
SAMPLESHEET = 'SampleSheet.csv'
BCL2FASTQ_REPORT = 'laneBarcode.html'


def getDatetime():
    try:
        d = datetime.datetime.now()
        return "{0}{1}{2}_{3}{4}{5}".format(
            d.day,
            d.month,
            d.year,
            d.hour,
            d.minute,
            d.second)
    except Exception as e:
        raise e


def getLogfile():
    try:
        d = getDatetime()
        logfile = os.open(os.path.join(
            WORKING_DIR, 'logfile-%s.log' % d), os.O_WRONLY | os.O_CREAT, 0o600)
        return logfile
    except Exception as e:
        raise e


def get_status_folder(file_status):
    if(not os.path.exists(file_status)):
        return False
    fs = open(file_status, 'r')
    status = fs.readline().strip()
    fs.close()

    return status


def get_run_details(args):
    try:
        if(os.path.exists(
                os.path.join(WORKING_DIR, args.runPath, SAMPLESHEET))):
            csv_file = open(os.path.join(WORKING_DIR, args.runPath, SAMPLESHEET), 'rb')
            ssheet = csv.reader(csv_file, delimiter=',')
            lines = OrderedDict([])
            key = ''
            not_null = [row for row in ssheet if len(row) > 0]
            for row in not_null:
                if(row[0].startswith('[')):
                    key = row[0].upper()
                    lines[key] = []
                else:
                    v = lines.get(key)
                    v.append(row)
                    lines[key] = v

            return lines
    except Exception as e:
        raise e


def get_bcl2fastq_report(args, fastq_path):
    try:

        if(os.path.exists(
                os.path.join(fastq_path, 'Reports'))):

            html = open(os.path.join(
                fastq_path,
                'Reports', 'html', 'index.html'), 'r').read()

            soup = BeautifulSoup(html, 'html.parser')
            for fr in soup.find_all('frame'):
                src = fr.get('src')
            src = src.replace('lane.html', BCL2FASTQ_REPORT)

            report = open(os.path.join(
                fastq_path,
                'Reports', 'html', src), 'r').read()

            soup = BeautifulSoup(report, 'html.parser')

            result = OrderedDict([])
            ncolums = 0

            hs = soup.find_all('h2')

            result['h2'] = [ele.text.strip() for ele in hs]

            tables = soup.find_all(id="ReportTable")
            for i, table in enumerate(tables):
                result['table-%i' % i] = OrderedDict([])
                for j, row in enumerate(table.find_all('tr')):
                    if('head' not in result['table-%i' % i]):
                        heads = row.find_all('th')
                        heads = [ele.text.strip() for ele in heads]
                        result['table-%i' % i]['head'] = heads
                        if(len(heads) > ncolums):
                            ncolums = len(heads)

                    cols = row.find_all('td')
                    cols = [ele.text.strip() for ele in cols]
                    if(len(cols) > 0):
                        result['table-%i' % i]['%i-col' % j] = cols
                        if(len(cols) > ncolums):
                            ncolums = len(cols)

            return ncolums, result
    except Exception as e:
        raise e


def rreplace(s, old, new, occurrence):

    li = s.rsplit(old, occurrence)
    return new.join(li)


def build_run_details_tex_table(args, data):
    if(data):
        tex_table = ''
        ncoluns = len(data['[DATA]'][0])

        # {|l|l|l|l|l|l|l|}
        columns_table = '{'
        for c in range(ncoluns):
            columns_table += '|l'
        columns_table += '|}'

        for key in data.keys():
            # HEADER
            values = data.get(key)
            tex_table += "\multicolumn{%s}{|c|}{%s} \\\\ \hline\n" % (
                ncoluns, key.replace('[', '').replace(']', ''))
            if(key == '[HEADER]'):
                for value in values:
                    tex_table += "%s & \multicolumn{%s}{l|}{%s} \\\\ \hline\n" % (
                        value[0].replace('_', '\_'), ncoluns - 1, value[1].replace('_', '\_'))
            # READS
            elif(key == '[READS]'):
                for value in values:
                    tex_table += "\multicolumn{%s}{|l|}{%s} \\\\ \hline\n" % (
                        ncoluns, value[0].replace('_', '\_'))
            # SETTINGS
            elif(key == '[SETTINGS]'):
                for value in values:
                    tex_table += "%s & \multicolumn{%s}{l|}{%s} \\\\ \hline\n" % (
                        value[0].replace('_', '\_'), ncoluns - 1, value[1].replace('_', '\_'))
            # DATA
            elif(key == '[DATA]'):
                for value in values:
                    tex_table += ''.join('%s & ' % v.replace('\\', '_').replace(
                        '_', '\_') for v in value)
                    tex_table = rreplace(tex_table, '&', ' ', 1)
                    tex_table += '\\\\ \hline\n'

        return columns_table, tex_table


def build_bcl2fastq_report_tex_table(args, fastq_path):
    ncoluns, data = get_bcl2fastq_report(args, fastq_path)

    if(data):
        tex_table = OrderedDict([])
        headers = data.get('h2')

        for i, head in enumerate(headers):

            if(head == 'Top Unknown Barcodes'):
                pass
            else:
                tb = data.get('table-%i' % i)

                tex = ''

                cols = len(tb['head'])

                tex += "\multicolumn{%s}{|c|}{%s} \\\\ \hline\n" % (cols, head)

                for key in tb.keys():
                    values = tb.get(key)

                    if(key == 'head'):
                        for v in values:
                            if(len(v.rsplit(" ", 1)) > 1):
                                v = "%s\\\\ %s" % (
                                    v.rsplit(" ", 1)[0], v.rsplit(" ", 1)[1])
                                line = "\\begin{tabular}[c]{@{}l@{}}%s\\end{tabular} &" % v.replace(
                                    '_', '\_').replace('%', '\%')
                                tex += line
                            else:
                                line = "%s &" % v
                                tex += line
                        tex = rreplace(tex, '&', ' ', 1)
                        tex += '\\\\ \hline\n'
                    else:
                        tex += ''.join('%s & ' % v.replace('_', '\_') for v in values)
                        tex = rreplace(tex, '&', ' ', 1)
                        tex += '\\\\ \hline\n'

            tex_table[head] = tex

        return tex_table


def check_analysed_folder(args, file_status):
    status = get_status_folder(file_status)
    if(status and status in ['emailed', 'running', 'completed']):
        return False
    if(not os.path.exists(file_status)):
        fs = open(file_status, 'w+')
        fs.write('checked\n')
        fs.close()

    return True


def run_blc2fastq(args, file_status, fastq_path, logfile):

    status = get_status_folder(file_status)

    if(status and status in ['converted']):
        return True

    if(os.path.exists(fastq_path)):
        return True

    cl = [
        '/usr/local/bin/bcl2fastq',
        '--runfolder-dir',
        args.runPath,
        '--output-dir',
        fastq_path]

    print('running blc2fastq')

    fs = open(file_status, 'w+')
    fs.write('running\n')
    fs.close()

    retProcess = subprocess.Popen(
        cl, 0, stdout=logfile, stderr=logfile, shell=False)
    retCode = retProcess.wait()
    if(retCode != 0):
        fs = open(file_status, 'w+')
        fs.write('error\n')
        fs.close()
        print(os.system('tail %s' % logfile))
        return False

    fs = open(file_status, 'w+')
    fs.write('converted\n')
    fs.close()

    print('finished')

    return True


def rename_fastq_file(args, fastq_path):
    try:
        fastq_files = {}

        for read in range(1, 3):
            npattern = 'L{0}_L00{0}_R{1}_{2}.'

            regex = '.*L00%d\\_R%d.*\\.gz\\Z(?ms)'
            lane = 'L00%d'
            if(args.sequencerName.upper() == 'NEXTSEQ'):
                lanes = 4
            else:
                lanes = 1
            for l in range(1, lanes + 1):  # NextSeq has 4 lanes
                clane = lane % l

                reobj = re.compile(regex % (l, read))

                files = [f for f in os.listdir(fastq_path) if reobj.match(f)]

                filedirs = [f for f in os.listdir(fastq_path) if os.path.isdir(
                    os.path.join(fastq_path, f))]

                for d in filedirs:
                    filelist = [f for f in os.listdir(
                        os.path.join(fastq_path, d)) if reobj.match(f)]
                    for f in filelist:
                        files.append(os.path.join(d, f))

                nfiles = []

                for i, f in enumerate(files):
                    name, ext = f.split('.', 1)

                    digits = len(str(i + 1))
                    if(digits == 1):
                        group = '00' + str(i + 1)
                    elif(digits == 2):
                        group = '0' + str(i + 1)
                    else:
                        group = str(i + 1)

                    f_npattern = npattern.format(l, read, group)
                    nname = f_npattern + ext

                    if(not os.path.islink(
                        os.path.join(
                            fastq_path, nname))):
                        os.symlink(
                            os.path.join(fastq_path, f),
                            os.path.join(fastq_path, nname))
                        nfiles.append(nname)

                if(nfiles):
                    if(clane in fastq_files):
                        values = fastq_files.get(clane)
                        for f in nfiles:
                            values.append(f)
                        fastq_files[clane] = values
                    else:
                        fastq_files[clane] = (nfiles)

        return fastq_files

    except Exception as e:
        raise e


def run_fastqc(args, file_status, fastq_path, logfile):
    status = get_status_folder(file_status)
    if(status and status in ['reported']):
        return True

    if(not os.path.exists(fastq_path)):
        return False

    fasta_files = rename_fastq_file(args, fastq_path)

    if(not fasta_files):
        return False

    # Check if there already is a fastq folder
    regex = '.*\\.html\\Z(?ms)'
    reobj = re.compile(regex)

    paths = [f for f in os.listdir(fastq_path) if reobj.match(f)]
    if(paths):
        return True

    if(args.sequencerName.upper() == 'NEXTSEQ'):
        lanes = 4
    else:
        lanes = 1

    for l in range(1, lanes + 1):
        lane = 'L00%d' % l
        files = fasta_files[lane]

        fileList = [' '.join(
            os.path.join(fastq_path, f) for f in files)]

        cl = ['/data/runs/FastQC/FastQC/fastqc --extract --casava -t 8 ' + fileList[0]]

        print('running fastqc')

        fs = open(file_status, 'w+')
        fs.write('running\n')
        fs.close()

        retProcess = subprocess.Popen(
            cl, 0, stdout=logfile, stderr=logfile, shell=True)
        retCode = retProcess.wait()
        if(retCode != 0):
            fs = open(file_status, 'w+')
            fs.write('error\n')
            fs.close()
            return False

        files = fasta_files[lane]
        for f in files:
            if(os.path.islink(os.path.join(fastq_path, f))):
                try:
                    os.unlink(os.path.join(fastq_path, f))
                except OSError as e:
                    'It was not possible to unlink the file \n%s. Error: %s' % (
                        os.path.join(fastq_path, f), e)

    fs = open(file_status, 'w+')
    fs.write('reported\n')
    fs.close()

    print('finished')

    return True


def compile_tex(args, file_status, fastq_path, logfile):
    status = get_status_folder(file_status)
    if(status and status in ['compiled']):
        return True

    images_dir = []
    reports_dir = []

    regex = '.*\\.html\\Z(?ms)'
    reobj = re.compile(regex)

    paths = [f for f in os.listdir(fastq_path) if reobj.match(f)]

    for path in paths:
        path_fastqc = path.split('.', 1)[0]

        s_image = os.path.join(fastq_path, path_fastqc, 'Images')

        report_dir = os.path.join(
            WORKING_DIR, args.runPath, REPORTS_PATH, path_fastqc)

        if(not os.path.exists(os.path.join(fastq_path, path_fastqc))):
            return False
        images_dir.append(s_image)
        reports_dir.append(report_dir)

    tex = open(os.path.join(WORKING_DIR, REPORT_FILE), 'r')
    rel = tex.read()
    tex.close()

    if(os.path.exists(os.path.join(WORKING_DIR, args.runPath, REPORTS_PATH))):
        shutil.rmtree(os.path.join(WORKING_DIR, args.runPath, REPORTS_PATH))

    os.mkdir(os.path.join(WORKING_DIR, args.runPath, REPORTS_PATH))

    data = get_run_details(args)
    tex_columns_table, tex_table_run_details = build_run_details_tex_table(args, data)

    tex_table_bcl2fastq_report = build_bcl2fastq_report_tex_table(args, fastq_path)

    for image_dir, report_dir in zip(images_dir, reports_dir):
        new_rel = rel.replace("$PATH$", image_dir)
        new_rel = new_rel.replace("$EQUIPAMENTO$", args.sequencerName)
        new_rel = new_rel.replace("$TABLECOLUMNS$", tex_columns_table)
        new_rel = new_rel.replace("$TABLECONTENTS$", tex_table_run_details)
        lane = report_dir.rsplit('_', 3)[1][-1]
        new_rel = new_rel.replace("$LANE$", lane)
        read = report_dir.rsplit('_', 2)[1]  # R1 or R2
        new_rel = new_rel.replace("$READ$", read)

        for i, key in enumerate(tex_table_bcl2fastq_report.keys()):
            char = chr(i + ord('A'))
            tex = tex_table_bcl2fastq_report.get(key)
            new_rel = new_rel.replace("$TABLE%sHEADER$" % char, key.encode('utf-8'))
            new_rel = new_rel.replace("$TABLE%sCONTENTS$" % char, tex.encode('utf-8'))

        os.mkdir(os.path.join(WORKING_DIR, args.runPath, REPORTS_PATH, report_dir))
        tex = open(
            os.path.join(WORKING_DIR, args.runPath, REPORTS_PATH, report_dir, REPORT_FILE), 'w+')
        tex.write(new_rel)
        tex.close()

        filename = '{0}-L00{1}-{2}'.format(REPORT_FILE.rsplit('.', 1)[0], lane, read)

        cl = [
            'pdflatex',
            '-output-directory',
            os.path.join(WORKING_DIR, args.runPath, REPORTS_PATH, report_dir),
            '--jobname=%s' % filename,
            os.path.join(WORKING_DIR, args.runPath, REPORTS_PATH, report_dir, REPORT_FILE)
        ]

        print('compiling tex')

        fs = open(file_status, 'w+')
        fs.write('running\n')
        fs.close()

        retProcess = subprocess.Popen(
            cl, 0, stdout=logfile, stderr=logfile, shell=False)
        retCode = retProcess.wait()
        if(retCode != 0):
            fs = open(file_status, 'w+')
            fs.write('error\n')
            fs.close()
            return False

    fs = open(file_status, 'w+')
    fs.write('compiled\n')
    fs.close()

    print('tex compiled')

    return True


def send_email():
    pass


def main():

    parser = argparse.ArgumentParser(description='Generate a PDF report with FastQC analysis')

    parser.add_argument(
        '--runPath', '-p', required=True,
        default=None, help='Path with the files of the run (default: %(default)s)')
    parser.add_argument(
        '--sequencerName', '-s',
        default='miseq',
        choices=['miseq', 'nextseq'],
        required=True,
        help='Sequencer name (default: %(default)s)')
    parser.add_argument(
        '--runName', '-r',
        default=None, help='Name of the run (default: %(default)s)')

    args = parser.parse_args()

    args.sequencerName = args.sequencerName.upper()

    file_status = os.path.join(WORKING_DIR, args.runPath, STATUS_FILE)

    if(not args.runName):
        if(args.runPath.endswith('/')):
            args.runName = os.path.join(WORKING_DIR, args.runPath).rsplit('/', 2)[-2]
        elif('/' in args.runPath):
            args.runName = os.path.join(WORKING_DIR, args.runPath).rsplit('/', 1)[-1]
        else:
            args.runName = os.path.join(WORKING_DIR, args.runPath)

    if(not os.path.exists(os.path.join(WORKING_DIR, args.runPath))):
        raise Exception(
            "Path of the run not found. \n %s" % os.path.join(WORKING_DIR, args.runPath))

    print('path exist')

    if(not check_analysed_folder(args, file_status)):
        raise Exception(
            'The folder has the status "%s". Execution aborted.' %
            get_status_folder(file_status).strip())

    print('path checked')

    fastq_path = ''

    logfile = getLogfile()

    fastq_path = os.path.join(WORKING_DIR, args.runPath, '%s_fastq/' % args.runName)
    if(not run_blc2fastq(args, file_status, fastq_path, logfile)):
        raise Exception("Error on bcl2fastq. Execution aborted.")

    print('converted')

    if(not run_fastqc(args, file_status, fastq_path, logfile)):
        raise Exception("Error on fastqc. Execution aborted.")

    print('reported')

    if(not compile_tex(args, file_status, fastq_path, logfile)):
        raise Exception("Error on compile tex. Execution aborted.")

    print('generated pdf')

    build_bcl2fastq_report_tex_table(args, fastq_path)

    try:
        os.remove(file_status)
    except Exception as e:
        raise e


if __name__ == '__main__':
    main()
