# -*- coding: utf-8 -*-
import time
import os
import uuid
import subprocess
import sys

from KBaseReport.KBaseReportClient import KBaseReport

from kb_kaiju.Utils.DataStagingUtils import DataStagingUtils
from kb_kaiju.Utils.OutputBuilder import OutputBuilder


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))
    sys.stdout.flush()


class KaijuUtil:

    def __init__(self, config, ctx):
        self.config = config
        self.ctx = ctx
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.threads = config['threads']
        self.SE_flag = 'SE'
        self.PE_flag = 'PE'
        

    def run_kaiju_and_krona(self, params):
        '''
        Main entry point for running kaiju + krona as a KBase App
        '''

        # 0) validate basic parameters
        method = 'run_kaiju_and_krona'
        required_params = ['workspace_name',
                           'input_refs',
                           'output_biom_name',
                           'tax_levels',
                           'db_type',
                           'seg_filter',
                           'min_match_length',
                           'greedy_run_mode',
                           'filter_percent',
                           'filter_unclassified',
                           'full_tax_path'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"' for method: '"+str(method)+"()'")

        if 'greedy_run_mode' in params and int(params['greedy_run_mode']) == 1:
            greedy_required_params = ['greedy_allowed_mismatches',
                                      'greedy_min_match_score',
                                      'greedy_max_e_value'
                                  ]
            for arg in greedy_required_params:
                if arg not in params or params[arg] == None or params[arg] == '':
                    raise ValueError ("Must define GREEDY MODE required param: '"+arg+"' for method: '"+str(method)+"()'")


        # 1) stage input data
        dsu = DataStagingUtils(self.config, self.ctx)
        staged_input = dsu.stage_input(params['input_refs'], 'fastq')
        input_dir = staged_input['input_dir']
        suffix = staged_input['folder_suffix']
        expanded_input = staged_input['expanded_input']

        output_dir = os.path.join(self.scratch, 'output_' + suffix)
        #plots_dir = os.path.join(self.scratch, 'plot_' + suffix)
        html_dir = os.path.join(self.scratch, 'html_' + suffix)

        log('Staged input directory: ' + input_dir)


        # 2) run Kaiju and make summaries
        kaiju_options = {'input_reads':               expanded_input,
                         'out_folder':                output_dir,
                         'tax_levels':                params['tax_levels'],
                         'db_type':                   params['db_type'],
                         'seg_filter':                params['seg_filter'],
                         'min_match_length':          params['min_match_length'],
                         'greedy_run_mode':           params['greedy_run_mode'],
                         'greedy_allowed_mismatches': params['greedy_allowed_mismatches'],
                         'greedy_min_match_score':    params['greedy_min_match_score'],
                         'filter_percent':            params['filter_percent'],
                         'filter_unclassified':       params['filter_unclassified'],
                         'full_tax_path':             params['full_tax_path'],
                         'threads':                   self.threads
                        }
        self.run_kaiju_batch (kaiju_options)

# HERE


        # 3) create Krona plots
        #self.build_checkM_lineage_wf_plots(input_dir, output_dir, plots_dir, all_seq_fasta_file, tetra_file)


        # 4) create Summary abundance plots 


        # 5) Package results
        #outputBuilder = OutputBuilder(output_dir, plots_dir, self.scratch, self.callback_url)
        #output_packages = self._build_output_packages(params, outputBuilder, input_dir)


        # 6) build the HTML report
        #os.makedirs(html_dir)
        #outputBuilder.build_html_output_for_lineage_wf(html_dir, params['input_ref'])
        #html_zipped = outputBuilder.package_folder(html_dir, 'report.html', 'Summarized report from CheckM')


        # 7) save report
        report_params = {'message': '',
                         #'direct_html_link_index': 0,
                         #'html_links': [html_zipped],
                         #'file_links': output_packages,
                         'report_object_name': 'kb_kaiju_report_' + str(uuid.uuid4()),
                         'workspace_name': params['workspace_name']
                         }

        kr = KBaseReport(self.callback_url)
        report_output = kr.create_extended_report(report_params)

        return {'report_name': report_output['name'],
                'report_ref':  report_output['ref']}


    def run_kaiju_batch(self, options, dropOutput=False):
        kaiju_options = {'input_reads':               expanded_input,
                         'out_folder':                output_dir,
                         'tax_levels':                params['tax_levels'],
                         'db_type':                   params['db_type'],
                         'seg_filter':                params['seg_filter'],
                         'min_match_length':          params['min_match_length'],
                         'greedy_run_mode':           params['greedy_run_mode'],
                         'greedy_allowed_mismatches': params['greedy_allowed_mismatches'],
                         'greedy_min_match_score':    params['greedy_min_match_score'],
                         'filter_percent':            params['filter_percent'],
                         'filter_unclassified':       params['filter_unclassified'],
                         'full_tax_path':             params['full_tax_path'],
                         'threads':                   self.threads
                        }


        input_reads = options['input_reads']
        for input_reads_item in input_reads:
            single_kaiju_run_options = kaiju_options
            single_kaiju_run_options['input_item'] = input_reads_item
            
            command = self._build_kaiju_command(options)
            log('Running: ' + ' '.join(command))

            log_output_file = None
            if dropOutput:  # if output is too chatty for STDOUT
                log_output_file = open(os.path.join(self.scratch, input_reads_item['name'] + '.out'), 'w')
                p = subprocess.Popen(command, cwd=self.scratch, shell=False, stdout=log_output_file, stderr=subprocess.STDOUT)
            else:
                p = subprocess.Popen(command, cwd=self.scratch, shell=False)
            exitCode = p.wait()

            if log_output_file:
                log_output_file.close()

            if (exitCode == 0):
                log('Executed command: ' + ' '.join(command) + '\n' +
                    'Exit Code: ' + str(exitCode))
            else:
                raise ValueError('Error running command: ' + ' '.join(command) + '\n' +
                                 'Exit Code: ' + str(exitCode))


    def _validate_kaiju_options(self, options):
        # 1st order required
        func_name = 'kaiju'
        required_opts = [ 'db_type',
                          'min_match_length',
                          'greedy_run_mode'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # 2nd order required
        if 'greedy_run_mode' in options and int(options['greedy_run_mode']) == 1:
            opt = 'greedy_allowed_mismatches'
            if opt not in options or int(options[opt]) < 1:
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()' if running in greedy_run_mode")

        # input file validation
        if not os.path.getsize(options['input_reads']['fwd_file']) > 0:
            raise ValueError ('missing or empty fwd reads file: '+options['input_reads']['fwd_file'])
        if options['input_reads']['type'] == self.PE_flag:
            if not os.path.getsize(options['input_reads']['rev_file']) > 0:
                raise ValueError ('missing or empty rev reads file: '+options['input_reads']['rev_file'])

        # db validation
        DB = 'KAIJU_DB_PATH'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaiju_options(self, command_list, options):
        kaiju_options = {'input_reads':               expanded_input,
                         'out_folder':                output_dir,
                         'tax_levels':                params['tax_levels'],
                         'db_type':                   params['db_type'],
                         'seg_filter':                params['seg_filter'],
                         'min_match_length':          params['min_match_length'],
                         'greedy_run_mode':           params['greedy_run_mode'],
                         'greedy_allowed_mismatches': params['greedy_allowed_mismatches'],
                         'greedy_min_match_score':    params['greedy_min_match_score'],
                         'filter_percent':            params['filter_percent'],
                         'filter_unclassified':       params['filter_unclassified'],
                         'full_tax_path':             params['full_tax_path'],
                         'threads':                   self.threads
                        }

        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_PATH'):
            command_list.append('-f')
            command_list.append(str(options.get('KAIJU_DB_PATH')))
        if options['input_reads'].get('fwd_file'):
            command_list.append('-i')
            command_list.append(str(options['input_reads'].get('fwd_file')))
        if options['input_reads'].get('type') == self.PE_flag:
            command_list.append('-j')
            command_list.append(str(options['input_reads'].get('rev_file')))
        if options.get('out_folder'):
            out_file = options['input_reads']['name']+'.kaiju'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)
        if int(options.get('seg_filter')) == 1:
            command_list.append('-x')
        if options.get('min_match_length'):
            command_list.append('-m')
            command_list.append(str(options.get('min_match_length')))
        if int(options.get('greedy_mode')) == 1:
            command_list.append('-a')
            command_list.append('greedy')
            if options.get('greedy_allowed_mismatches'):
                command_list.append('-e')
                command_list.append(str(options.get('greedy_allowed_mismatches')))
            if options.get('greedy_min_match_score'):
                command_list.append('-s')
                command_list.append(str(options.get('greedy_min_match_score')))

        if options.get('threads'):
            command_list.append('-z')
            command_list.append(str(options.get('threads')))
        if options.get('verbose'):
            command_list.append('-v')
        

    def _build_kaiju_command(self, options, verbose=True):

        KAIJU_BIN_DIR  = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_BIN      = os.path.join(KAIJU_BIN_DIR, 'kaiju')
        KAIJU_DB_DIR   = os.path.join(os.path.sep, 'data', 'kb_kaiju', 'kaijudb', options['db_type'])
        options['verbose'] = verbose
        if self.threads and self.threads > 1:
            options['threads'] = self.threads

        options['KAIJU_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        #options['KAIJU_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')  # don't need for kaiju cmd

        if options['db_type'] == 'kaiju_index':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db.fmi')
        elif options['db_type'] == 'kaiju_index_pg':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db.fmi')
        elif options['db_type'] == 'kaiju_index_nr':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr.fmi')
        elif options['db_type'] == 'kaiju_index_nr_euk':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr_euk.fmi')
        else:
            raise ValueError ('bad db_type: '+options['db_type']+' (must be one of "kaiju_index", "kaiju_index_pg", "kaiju_index_nr", "kaiju_index_nr_euk")')

        self._validate_kaiju_options(options)

        command = [KAIJU_BIN]
        self._process_kaiju_options(command, options)

        return command


    def _validate_kaijuReport_options(self, options,
                          checkBin=False,
                          checkOut=False,
                          checkPlots=False,
                          checkTetraFile=False,
                          subcommand=''):
        # Note: we can, maybe should, add additional checks on the contents of the folders here
        if checkBin and 'bin_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without bin_folder option set')
        if checkOut and 'out_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without bin_folder option set')
        if checkPlots and 'plots_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without plots_folder option set')
        if checkTetraFile and 'tetra_file' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without tetra_file option set')


    def _process_kaijuReport_options(self, command_list, options):
        if options.get('thread'):
            command_list.append('-t')
            command_list.append(str(options.get('thread')))

        if options.get('quiet') and str(options.get('quiet')) == '1':
            command_list.append('--quiet')


    def _build_kaijuReport_command(self, options):

        KAIJU_BIN_DIR    = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_REPORT_BIN = os.path.join(KAIJU_BIN_DIR, 'kaijuReport')

        self._validate_kaijuReport_options(options)

        command = [KAIJU_REPORT_BIN]
        self._process_kaijuReport_options(command, options)

        return command


    def _validate_kaiju2krona_options(self, options,
                          checkBin=False,
                          checkOut=False,
                          checkPlots=False,
                          checkTetraFile=False,
                          subcommand=''):
        # Note: we can, maybe should, add additional checks on the contents of the folders here
        if checkBin and 'bin_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without bin_folder option set')
        if checkOut and 'out_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without bin_folder option set')
        if checkPlots and 'plots_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without plots_folder option set')
        if checkTetraFile and 'tetra_file' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without tetra_file option set')


    def _process_kaiju2krona_options(self, command_list, options):
        if options.get('thread'):
            command_list.append('-t')
            command_list.append(str(options.get('thread')))

        if options.get('quiet') and str(options.get('quiet')) == '1':
            command_list.append('--quiet')


    def _build_kaiju2krona_command(self, options):

        KAIJU_BIN_DIR   = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU2KRONA_BIN = os.path.join(KAIJU_BIN_DIR, 'kaiju2krona')

        self._validate_kaiju2krona_options(options)

        command = [KAIJU2KRONA_BIN]
        self._process_kaiju2krona_options(command, options)

        return command


    def _validate_kronaImport_options(self, options,
                          checkBin=False,
                          checkOut=False,
                          checkPlots=False,
                          checkTetraFile=False,
                          subcommand=''):
        # Note: we can, maybe should, add additional checks on the contents of the folders here
        if checkBin and 'bin_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without bin_folder option set')
        if checkOut and 'out_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without bin_folder option set')
        if checkPlots and 'plots_folder' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without plots_folder option set')
        if checkTetraFile and 'tetra_file' not in options:
            raise ValueError('cannot run checkm ' + subcommand + ' without tetra_file option set')


    def _process_kronaImport_options(self, command_list, options):
        if options.get('thread'):
            command_list.append('-t')
            command_list.append(str(options.get('thread')))

        if options.get('quiet') and str(options.get('quiet')) == '1':
            command_list.append('--quiet')


    def _build_kronaImport_command(self, options):

        KRONA_BIN_DIR    = os.path.join(os.path.sep, 'kb', 'module', 'Krona', 'bin')
        KRONA_IMPORT_BIN = os.path.join(KRONA_BIN_DIR, 'ktImportText')

        self._validate_kronaImport_options(options)

        command = [KRONA_IMPORT_BIN]
        self._process_kronaImport_options(command, options)

        return command


    def _build_output_packages(self, params, outputBuilder, input_dir):

        output_packages = []

        #if 'save_output_dir' in params and str(params['save_output_dir']) == '1':
        if True:
            log('packaging full output directory')
            zipped_output_file = outputBuilder.package_folder(outputBuilder.output_dir, 'full_output.zip',
                                                              'Full output of CheckM')
            output_packages.append(zipped_output_file)
        else:  # ADD LATER?
            log('not packaging full output directory, selecting specific files')
            crit_out_dir = os.path.join(self.scratch, 'critical_output_' + os.path.basename(input_dir))
            os.makedirs(crit_out_dir)
            zipped_output_file = outputBuilder.package_folder(outputBuilder.output_dir, 'selected_output.zip',
                                                              'Selected output from the CheckM analysis')
            output_packages.append(zipped_output_file)


        if 'save_plots_dir' in params and str(params['save_plots_dir']) == '1':
            log('packaging output plots directory')
            zipped_output_file = outputBuilder.package_folder(outputBuilder.plots_dir, 'plots.zip',
                                                              'Output plots from CheckM')
            output_packages.append(zipped_output_file)
        else:
            log('not packaging output plots directory')

        return output_packages
