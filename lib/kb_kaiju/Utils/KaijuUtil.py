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
        self.workspace_url = config['workspace-url']
        self.scratch = config['scratch']
        self.threads = config['threads']
        self.suffix = str(int(time.time() * 1000))
        self.SE_flag = 'SE'
        self.PE_flag = 'PE'
        self.dsu_client = DataStagingUtils(self.config, self.ctx)

        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)


    def run_kaiju_with_krona(self, params):
        '''
        Main entry point for running kaiju + krona as a KBase App
        '''

        # 0) validate basic parameters and set defaults
        params = self.validate_run_kaiju_with_krona_params(params)

        # 1) expand input members that are sets
        expanded_input = self.dsu_client.expand_input(params['input_refs'])


        # 2) establish output folders
        output_dir = os.path.join(self.scratch, 'output_' + str(self.suffix))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        html_dir = os.path.join(self.scratch, 'html_' + str(self.suffix))
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)

        kaiju_output_folder = os.path.join(output_dir, 'kaiju_output')
        if not os.path.exists(kaiju_output_folder):
            os.makedirs(kaiju_output_folder)

        kaijuReport_output_folder = os.path.join(output_dir, 'kaiju_report')
        if not os.path.exists(kaijuReport_output_folder):
            os.makedirs(kaijuReport_output_folder)

        kaijuReport_PerSamplePlots_output_folder = os.path.join(output_dir, 'kaiju_report_per_sample_plots')
        kaijuReport_StackedBarPlots_output_folder = os.path.join(output_dir, 'kaiju_report_stacked_bar_plots')
        kaijuReport_StackedAreaPlots_output_folder = os.path.join(output_dir, 'kaiju_report_stacked_area_plots')
        if not os.path.exists(kaijuReport_PerSamplePlots_output_folder):
            os.makedirs(kaijuReport_PerSamplePlots_output_folder)
        if not os.path.exists(kaijuReport_StackedBarPlots_output_folder):
            os.makedirs(kaijuReport_StackedBarPlots_output_folder)
        if not os.path.exists(kaijuReport_StackedAreaPlots_output_folder):
            os.makedirs(kaijuReport_StackedAreaPlots_output_folder)

        krona_output_folder = os.path.join(output_dir, 'krona_data')
        if not os.path.exists(krona_output_folder):
            os.makedirs(krona_output_folder)


        # 3) instantiate OutputBuilder
        output_folders = [ { 'name': 'kaiju_classifications',
                             'desc': 'Kaiju Classification',
                             'path': kaiju_output_folder
                         },
                           { 'name': 'kaiju_summaries',
                             'desc': 'Kaiju Summaries',
                             'path': kaijuReport_output_folder
                         },
                           { 'name': 'krona_data',
                             'desc': 'Krona Data',
                             'path': krona_output_folder
                         },
                           { 'name': 'stacked_bar_abundance_plots_PNG+PDF',
                             'desc': 'Stacked Bar Abundance Plots (PNG + PDF)',
                             'path': kaijuReport_StackedBarPlots_output_folder
                         }
                           #{ 'name': 'per_sample_ranked_abundance_plots_PNG+PDF',
                           #  'desc': 'Per Sample Ranked Abundance Plots (PNG + PDF)',
                           #  'path': kaijuReport_PerSamplePlots_output_folder
                           #}
                      ]
        build_area_plots_flag = False
        # DEBUG
        #if len(expanded_input) > 1:
        #    build_area_plots_flag = True
        #    output_folders.append ({ 'name': 'stacked_area_abundance_plots_PNG+PDF',
        #                             'desc': 'Stacked Area Abundance Plots (PNG + PDF)',
        #                             'path': kaijuReport_StackedAreaPlots_output_folder
        #                           })
        self.outputBuilder_client = OutputBuilder(output_folders, self.scratch, self.callback_url, self.workspace_url)


        # 4) run Kaiju in batch (download happens one-by-one and then deleted to save space)
        kaiju_options = {'input_reads':               expanded_input,
                         'out_folder':                kaiju_output_folder,
                         'subsample_percent':         params['subsample_percent'],
                         'subsample_replicates':      params['subsample_replicates'],
                         'subsample_seed':            params['subsample_seed'],
                         'tax_levels':                params['tax_levels'],
                         'db_type':                   params['db_type'],
                         'seg_filter':                params['seg_filter'],
                         'min_match_length':          params['min_match_length'],
                         'greedy_run_mode':           params['greedy_run_mode'],
                         'greedy_allowed_mismatches': params['greedy_allowed_mismatches'],
                         'greedy_min_match_score':    params['greedy_min_match_score'],
                         'threads':                   self.threads
                        }
        expanded_input = self.run_kaiju_batch (kaiju_options)  # revise expanded input with subsamples


        # 5) create Summary Reports in batch
        kaijuReport_options = {'input_reads':               expanded_input,
                               'in_folder':                 kaiju_output_folder,
                               'out_folder':                kaijuReport_output_folder,
                               'tax_levels':                params['tax_levels'],
                               'db_type':                   params['db_type'],
                               'filter_percent':            params['filter_percent'],
                               'filter_unclassified':       params['filter_unclassified'],
                               'full_tax_path':             params['full_tax_path']
                           }
        self.run_kaijuReport_batch (kaijuReport_options)


        # 6) create Summary Report plots in batch
        kaijuReportPlots_options = {'input_reads':                   expanded_input,
                                    'in_folder':                     kaijuReport_output_folder,
                                    'stacked_bar_plots_out_folder':  kaijuReport_StackedBarPlots_output_folder,
                                    #'per_sample_plots_out_folder':   kaijuReport_PerSamplePlots_output_folder,
                                    'tax_levels':                    params['tax_levels'],
                                    'sort_taxa_by':                  params['sort_taxa_by']
                                    #'filter_percent':            params['filter_percent'],
                                    #'filter_unclassified':       params['filter_unclassified'],
                                    #'full_tax_path':             params['full_tax_path']
                                }
        if build_area_plots_flag:
            kaijuReportPlots_options['stacked_area_plots_out_folder'] = kaijuReport_StackedAreaPlots_output_folder
        kaijuReport_plot_files = self.run_kaijuReportPlots_batch (kaijuReportPlots_options)


        # 7) create HTML Summary Reports in batch
        kaijuReportPlotsHTML_options = {'input_reads':             expanded_input,
                                        'summary_folder':          kaijuReport_output_folder,
                                        'stacked_bar_plot_files':  kaijuReport_plot_files['stacked_bar_plot_files'],
                                        #'per_sample_plot_files':   kaijuReport_plot_files['per_sample_plot_files'],
                                        'out_folder':              html_dir,
                                        'tax_levels':              params['tax_levels']
        }
        if build_area_plots_flag:
            kaijuReportPlotsHTML_options['stacked_area_plot_files'] = kaijuReport_plot_files['stacked_area_plot_files']
        html_plot_pages = self.run_kaijuReportPlotsHTML_batch (kaijuReportPlotsHTML_options)


        # 8) create Krona plots
        krona_options = {'input_reads':               expanded_input,
                         'in_folder':                 kaiju_output_folder,
                         'out_folder':                krona_output_folder,
                         'html_folder':               html_dir,
                         'db_type':                   params['db_type']
                     }
        html_krona_pages = self.run_krona_batch (krona_options)


        # 9) Package results
        output_packages = self._build_output_packages(params, self.outputBuilder_client)


        # 10) add top nav to html pages and build the HTML report
        html_pages = []
        html_pages.extend(html_plot_pages['bar'])
        if build_area_plots_flag:
            html_pages.extend(html_plot_pages['area'])
        #html_pages.extend(html_plot_pages['per_sample'])
        html_pages.extend(html_krona_pages)
        self.outputBuilder_client.add_top_nav(html_pages)
        #report_html_file = 'kaiju_html_plots.zip'  # fails
        #report_html_file = 'kaiju_plots.html'  # fails
        report_html_file = html_pages[0]['local_path']  # works
        report_html_desc = 'Kaiju abundance and Krona plots'
        html_zipped = self.outputBuilder_client.package_folder(html_dir, report_html_file, report_html_desc)


        """
        # 11) save biom output for each tax level
        generated_biom_objs = []
        timestamp_epoch = time.time()
        for tax_level in params['tax_levels']:
            obj_name = params['output_biom_name']
            if len(params['tax_levels']) != 1:
                obj_name += '-'+tax_level

            generate_biom_options = {'tax_level':       tax_level,
                                     'db_type':         params['db_type'],
                                     'input_reads':     expanded_input,
                                     'in_folder':       kaiju_output_folder,
                                     'workspace_name':  params['workspace_name'],
                                     'output_obj_name': obj_name,
                                     'timestamp_epoch': timestamp_epoch
                                 }
            biom_obj_ref = self.outputBuilder_client.generate_sparse_biom1_0_matrix(self.ctx, generate_biom_options)
            generated_biom_objs.append({'ref': biom_obj_ref,
                                        'description': 'Kaiju Taxonomic Classification at '+tax_level+' Level.  BIOM format'})
        """


        # 12) save report
        report_params = {'message': '',
                         #'objects_created': generated_biom_objs,
                         'objects_created': [],
                         'direct_html_link_index': 0,
                         'html_links': [html_zipped],
                         'file_links': output_packages,
                         'report_object_name': 'kb_kaiju_report_' + str(uuid.uuid4()),
                         'workspace_name': params['workspace_name']
                         }

        kr = KBaseReport(self.callback_url)
        report_output = kr.create_extended_report(report_params)

        returnVal = {'report_name': report_output['name'],
                     'report_ref':  report_output['ref']}
        return returnVal


    def run_proc(self, command, log_output_file=None):
        log('Running: ' + ' '.join(command))

        if log_output_file:  # if output is too chatty for STDOUT
            open (log_output_file, 'w')
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
        return exitCode


    def validate_run_kaiju_with_krona_params(self, params):
        method = 'run_kaiju_with_krona'

        # base required params
        required_params = ['workspace_name',
                           'input_refs',
                           #'output_biom_name',
                           'tax_levels',
                           'db_type',
                           'filter_percent',
                           'seg_filter',
                           'min_match_length',
                           'greedy_run_mode',
                           'sort_taxa_by'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"' for method: '"+str(method)+"()'")

        # custom mode params
        if 'greedy_run_mode' in params and int(params['greedy_run_mode']) == 1:
            greedy_required_params = ['greedy_allowed_mismatches',
                                      'greedy_min_match_score',
                                      'greedy_max_e_value'
                                  ]
            for arg in greedy_required_params:
                if arg not in params or params[arg] == None or params[arg] == '':
                    raise ValueError ("Must define GREEDY MODE required param: '"+arg+"' for method: '"+str(method)+"()'")

        # default vals for not required params
        default_param_vals = {'subsample_percent': 10,
                              'subsample_replicates': 1,
                              'subsample_seed': int(self.suffix),
                              'filter_unclassified': 1,
                              'full_tax_path': 0
                          }
        for arg in default_param_vals.keys():
            if arg not in params or params[arg] == None or params[arg] == '':
                params[arg] = default_param_vals[arg]

        # check math
        total_perc = float(params['subsample_percent']) * int(params['subsample_replicates'])
        if total_perc > 100:
            raise ValueError ("Subsample is non-overlapping, so too many subsample replicates "+str(params['subsample_replicates'])+" at subsample percent: "+str(params['subsample_perc'])+" (replicates * percent = "+str(total_perc)+" > 100)")

        # adjust param values by flag
        tax_levels_all = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        for tax_level in params['tax_levels']:
            if tax_level == 'ALL':
                params['tax_levels'] = tax_levels_all
                break
            elif tax_level not in tax_levels_all:
                raise ValueError ("Bad tax level "+tax_level)

        # make sure min and max vals not exceeded (input widget constraints not reliable)
        limit_vals = {'filter_percent':            {'min': 0, 'max': 10 },
                      'subsample_percent':         {'min': 1, 'max': 100 },
                      'subsample_replicates':      {'min': 1 },
                      'min_match_length':          {'min': 9 },
                      'greedy_allowed_mismatches': {'min': 1 },
                      'greedy_min_match_score':    {'min': 10},
                      'greedy_max_e_value':        {'min': 0.0, 'max': 1.0}
                  }
        bad_vals_msgs = []
        for arg in params.keys():
            if arg in limit_vals:
                if 'min' in limit_vals[arg] and float(params[arg]) < limit_vals[arg]['min']:
                    bad_vals_msgs.append('Value less than minimum for parameter '+arg+' ('+str(param[arg])+' < '+str(limit_vals[arg]['min']))
                elif 'max' in limit_vals[arg] and float(params[arg]) > limit_vals[arg]['max']:
                    bad_vals_msgs.append('Value greater than maximum for parameter '+arg+' ('+str(param[arg])+' > '+str(limit_vals[arg]['max']))
        if len(bad_vals_msgs) > 0:
            raise ValueError ("\n".join(bad_vals_msgs)+"\n")

        # return adjusted params
        return params


    def run_kaiju_batch(self, options, dropOutput=False):
        new_expanded_input = []

        input_reads = options['input_reads']
        for input_reads_item in input_reads:

            # download and subsample reads
            staged_input = self.dsu_client.stage_input(input_item =           input_reads_item,
                                                       subsample_percent =    int(options['subsample_percent']),
                                                       subsample_replicates = int(options['subsample_replicates']),
                                                       subsample_seed =       int(options['subsample_seed']),
                                                       fasta_file_extension = 'fastq')
            #input_dir = staged_input['input_dir']
            replicate_input = staged_input['replicate_input']
            new_expanded_input.extend(replicate_input)  # revise expanded input to replicates

            # run for each replicate
            for input_reads_item_replicate in replicate_input:
                single_kaiju_run_options = options
                single_kaiju_run_options['input_item'] = input_reads_item_replicate

                print ("REPLICATE: "+str(input_reads_item_replicate))

                log_output_file = None
                if dropOutput:  # if output is too chatty for STDOUT
                    log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kaiju' + '.stdout')

                command = self._build_kaiju_command(single_kaiju_run_options)
                self.run_proc (command, log_output_file)

                # remove input file to free up disk
                os.remove(input_reads_item_replicate['fwd_file'])
                if input_reads_item_replicate['type'] == self.PE_flag:
                    os.remove(input_reads_item_replicate['rev_file'])

        return new_expanded_input


    def run_kaijuReport_batch(self, options, dropOutput=False):
        input_reads = options['input_reads']
        for input_reads_item in input_reads:
            for tax_level in options['tax_levels']:
                single_kaijuReport_run_options = options
                single_kaijuReport_run_options['input_item'] = input_reads_item
                single_kaijuReport_run_options['tax_level'] = tax_level

                log_output_file = None
                if dropOutput:  # if output is too chatty for STDOUT
                    log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kaijuReport' + '.stdout')

                command = self._build_kaijuReport_command(single_kaijuReport_run_options)
                self.run_proc (command, log_output_file)


    def run_kaijuReportPlots_batch(self, options):
        input_reads = options['input_reads']
        per_sample_plot_files   = dict()
        stacked_bar_plot_files  = dict()
        stacked_area_plot_files = dict()

        for tax_level in options['tax_levels']:

            # per sample plots
            if 'per_sample_plots_out_folder' in options:
                per_sample_plot_files[tax_level] = dict()
                for input_reads_item in input_reads:
                    single_kaijuReportPlots_options = options
                    single_kaijuReportPlots_options['input_item'] = input_reads_item
                    single_kaijuReportPlots_options['tax_level'] = tax_level

                    per_sample_plot_files[tax_level][input_reads['name']] = self.outputBuilder_client.generate_kaijuReport_PerSamplePlots(single_kaijuReportPlots_options)

            # stacked bar plots
            if 'stacked_bar_plots_out_folder' in options:
                kaijuReportPlots_options = options
                kaijuReportPlots_options['stacked_plots_out_folder'] = options['stacked_bar_plots_out_folder']
                kaijuReportPlots_options['tax_level'] = tax_level
                kaijuReportPlots_options['plot_type'] = 'bar'
                stacked_bar_plot_files[tax_level] = self.outputBuilder_client.generate_kaijuReport_StackedPlots(kaijuReportPlots_options)

            # stacked area plots
            if 'stacked_area_plots_out_folder' in options:
                kaijuReportPlots_options = options
                kaijuReportPlots_options['stacked_plots_out_folder'] = options['stacked_area_plots_out_folder']
                kaijuReportPlots_options['tax_level'] = tax_level
                kaijuReportPlots_options['plot_type'] = 'area'
                stacked_area_plot_files[tax_level] = self.outputBuilder_client.generate_kaijuReport_StackedPlots(kaijuReportPlots_options)

        return {'per_sample_plot_files': per_sample_plot_files,
                'stacked_bar_plot_files': stacked_bar_plot_files,
                'stacked_area_plot_files': stacked_area_plot_files
            }


    def run_kaijuReportPlotsHTML_batch(self, options):
        out_html_folder = options['out_folder']
        out_html_files = dict()

        if 'stacked_bar_plot_files' in options:
            out_html_files['bar'] = self.outputBuilder_client.build_html_for_kaijuReport_StackedPlots(
                options['input_reads'],
                options['summary_folder'],
                out_html_folder,
                'bar',
                options['tax_levels'],
                options['stacked_bar_plot_files']
            )

        if 'stacked_area_plot_files' in options:
            out_html_files['area'] = self.outputBuilder_client.build_html_for_kaijuReport_StackedPlots(
                options['input_reads'],
                options['summary_folder'],
                out_html_folder,
                'area',
                options['tax_levels'],
                options['stacked_area_plot_files']
            )

        if 'per_sample_plot_files' in options:
            out_html_files['per_sample'] = self.outputBuilder_client.build_html_for_kaijuReport_StackedPlots(
                options['input_reads'],
                options['summary_folder'],
                out_html_folder,
                options['tax_levels'],
                options['per_sample_plot_files']
            )

        return out_html_files


    def run_krona_batch(self, options, dropOutput=False):
        out_html_files = []
        input_reads = options['input_reads']
        for input_reads_item in input_reads:

            # kaiju2krona
            single_kaiju2krona_run_options = options
            single_kaiju2krona_run_options['input_item'] = input_reads_item

            log_output_file = None
            if dropOutput:  # if output is too chatty for STDOUT
                log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kaiju2krona' + '.stdout')

            command = self._build_kaiju2krona_command(single_kaiju2krona_run_options)
            self.run_proc (command, log_output_file)

            # kronaImport
            single_kronaImport_run_options = options
            single_kronaImport_run_options['input_item'] = input_reads_item

            log_output_file = None
            if dropOutput:  # if output is too chatty for STDOUT
                log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kronaImport' + '.stdout')

            command = self._build_kronaImport_command(single_kronaImport_run_options)
            self.run_proc (command, log_output_file)

            # return file info
            local_html_path = input_reads_item['name']+'.krona.html'
            html_path = os.path.join (options['html_folder'], local_html_path)
            out_html_files.append({'type': 'krona',
                                   'name': input_reads_item['name']+' Krona',
                                   'local_path': local_html_path,
                                   'abs_path': html_path
                               })

        return out_html_files


    def _validate_kaiju_options(self, options):
        # 1st order required
        func_name = 'kaiju'
        required_opts = [ 'input_item',
                          'out_folder',
                          'db_type',
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
        if not os.path.getsize(options['input_item']['fwd_file']) > 0:
            raise ValueError ('missing or empty fwd reads file: '+options['input_item']['fwd_file'])
        if options['input_item']['type'] == self.PE_flag:
            if not os.path.getsize(options['input_item']['rev_file']) > 0:
                raise ValueError ('missing or empty rev reads file: '+options['input_item']['rev_file'])

        # db validation
        DB = 'KAIJU_DB_PATH'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaiju_options(self, command_list, options):
        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_PATH'):
            command_list.append('-f')
            command_list.append(str(options.get('KAIJU_DB_PATH')))
        if options['input_item'].get('fwd_file'):
            command_list.append('-i')
            command_list.append(str(options['input_item'].get('fwd_file')))
        if options['input_item'].get('type') == self.PE_flag:
            command_list.append('-j')
            command_list.append(str(options['input_item'].get('rev_file')))
        if options.get('out_folder'):
            out_file = options['input_item']['name']+'.kaiju'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)
        if int(options.get('seg_filter')) == 1:
            command_list.append('-x')
        if options.get('min_match_length'):
            command_list.append('-m')
            command_list.append(str(options.get('min_match_length')))
        if int(options.get('greedy_run_mode')) == 1:
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
        KAIJU_DB_DIR   = os.path.join(os.path.sep, 'data', 'kaijudb', options['db_type'])

        options['verbose'] = verbose
        if self.threads and self.threads > 1:
            options['threads'] = self.threads
        options['KAIJU_DB_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        #options['KAIJU_DB_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')  # don't need for kaiju cmd
        if options['db_type'] == 'refseq':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_refseq.fmi')
        elif options['db_type'] == 'progenomes':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_progenomes.fmi')
        elif options['db_type'] == 'nr':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr.fmi')
        elif options['db_type'] == 'nr_euk':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr_euk.fmi')
        else:
            raise ValueError ('bad db_type: '+options['db_type']+' (must be one of "refseq", "progenomes", "nr", "nr_euk")')

        self._validate_kaiju_options(options)
        command = [KAIJU_BIN]
        self._process_kaiju_options(command, options)
        return command


    def _validate_kaijuReport_options(self, options):
        # 1st order required
        func_name = 'kaiju2table'
        required_opts = [ 'in_folder',
                          'input_item',
                          'out_folder',
                          'db_type',
                          'tax_level'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # input file validation
        in_file = os.path.join(options['in_folder'], options['input_item']['name']+'.kaiju')
        if not os.path.getsize(in_file) > 0:
            raise ValueError ('missing or empty kaiju classification file: '+in_file)

        # db validation
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NAMES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaijuReport_options(self, command_list, options):
        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_NAMES'):
            command_list.append('-n')
            command_list.append(str(options.get('KAIJU_DB_NAMES')))
        if options.get('tax_level'):
            command_list.append('-r')
            command_list.append(str(options.get('tax_level')))
        if options.get('out_folder'):
            out_file = options['input_item']['name']+'-'+str(options.get('tax_level'))+'.kaijuReport'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)
        if options.get('filter_percent'):
            command_list.append('-m')
            command_list.append(str(options.get('filter_percent')))
        if int(options.get('filter_unclassified')) == 1:
            command_list.append('-u')
        if int(options.get('full_tax_path')) == 1:
            command_list.append('-p')
        if options.get('in_folder'): # needs to be last
            in_path = options['in_folder']
            command_list.append(in_path)

    def _build_kaijuReport_command(self, options):
        KAIJU_BIN_DIR    = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_REPORT_BIN = os.path.join(KAIJU_BIN_DIR, 'kaiju2table')
        KAIJU_DB_DIR     = os.path.join(os.path.sep, 'data', 'kaijudb', options['db_type'])

        options['KAIJU_DB_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        options['KAIJU_DB_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')

        self._validate_kaijuReport_options(options)
        command = [KAIJU_REPORT_BIN]
        self._process_kaijuReport_options(command, options)
        return command


    def _validate_kaiju2krona_options(self, options):
        # 1st order required
        func_name = 'kaiju2krona'
        required_opts = [ 'in_folder',
                          'input_item',
                          'out_folder',
                          'db_type'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # input file validation
        in_file = os.path.join(options['in_folder'], options['input_item']['name']+'.kaiju')
        if not os.path.getsize(in_file) > 0:
            raise ValueError ('missing or empty kaiju classification file: '+in_file)

        # db validation
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NAMES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaiju2krona_options(self, command_list, options):
        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_NAMES'):
            command_list.append('-n')
            command_list.append(str(options.get('KAIJU_DB_NAMES')))
        if options.get('in_folder'):
            in_file = options['input_item']['name']+'.kaiju'
            in_path = os.path.join(options['in_folder'], in_file)
            command_list.append('-i')
            command_list.append(in_path)
        if options.get('out_folder'):
            out_file = options['input_item']['name']+'.krona'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)


    def _build_kaiju2krona_command(self, options):
        KAIJU_BIN_DIR     = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_2_KRONA_BIN = os.path.join(KAIJU_BIN_DIR, 'kaiju2krona')
        KAIJU_DB_DIR      = os.path.join(os.path.sep, 'data', 'kaijudb', options['db_type'])

        options['KAIJU_DB_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        options['KAIJU_DB_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')

        self._validate_kaiju2krona_options(options)
        command = [KAIJU_2_KRONA_BIN]
        self._process_kaiju2krona_options(command, options)
        return command


    def _validate_kronaImport_options(self, options):
        # 1st order required
        func_name = 'kronaImport'
        required_opts = [ 'html_folder',
                          'input_item',
                          'out_folder',
                          'db_type'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # input file validation
        in_file = os.path.join(options['out_folder'], options['input_item']['name']+'.krona')
        if not os.path.getsize(in_file) > 0:
            raise ValueError ('missing or empty krona input file: '+in_file)


    def _process_kronaImport_options(self, command_list, options):
        if options.get('html_folder'):
            html_file = options['input_item']['name']+'.krona.html'
            html_path = os.path.join (str(options.get('html_folder')), html_file)
            command_list.append('-o')
            command_list.append(html_path)
        if options.get('out_folder'):
            in_file = options['input_item']['name']+'.krona'
            in_path = os.path.join(options['out_folder'], in_file)
            command_list.append(in_path)


    def _build_kronaImport_command(self, options):
        KRONA_BIN_DIR    = os.path.join(os.path.sep, 'usr', 'local', 'bin')
        KRONA_IMPORT_BIN = os.path.join(KRONA_BIN_DIR, 'ktImportText')
        #KRONA_BIN_DIR    = os.path.join(os.path.sep, 'kb', 'module', 'Krona', 'KronaTools', 'scripts')
        #KRONA_IMPORT_BIN = os.path.join(KRONA_BIN_DIR, 'ImportText.pl')

        self._validate_kronaImport_options(options)
        command = [KRONA_IMPORT_BIN]
        self._process_kronaImport_options(command, options)
        return command


    def _build_output_packages(self, params, outputBuilder):

        output_packages = []
        for output_folder in outputBuilder.output_folders:
            if 'skip_output_dirs' in params and output_folder['name'] in params['skip_output_dirs']:
                log('skipping output directory '+output_folder['name'])
                continue

            log('packaging output directory '+output_folder['name'])
            zipped_output_file = outputBuilder.package_folder(output_folder['path'],
                                                              output_folder['name']+'.zip',
                                                              output_folder['desc'])
            output_packages.append(zipped_output_file)

        return output_packages
