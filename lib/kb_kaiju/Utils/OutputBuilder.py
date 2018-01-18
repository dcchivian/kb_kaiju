import os
import shutil
import ast
import sys
import time

from DataFileUtil.DataFileUtilClient import DataFileUtil


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))
    sys.stdout.flush()


class OutputBuilder(object):
    '''
    Constructs the output HTML report and artifacts based on Kaiju and Krona
    runs.  This includes creating matplotlib graphs of the summaries and
    modifying the Krona HTML to offer tabbed href links between html pages
    '''

    def __init__(self, output_dirs, scratch_dir, callback_url):
        self.output_dirs = output_dirs
        self.scratch = scratch_dir
        self.callback_url = callback_url

        # leave out light colors
        self.no_light_color_names = [
            #'aliceblue',
            'aqua',
            'aquamarine',
            #'azure',
            #'beige',
            #'bisque',
            #'blanchedalmond',
            'blue',
            'blueviolet',
            'brown',
            'burlywood',
            'cadetblue',
            'chartreuse',
            'chocolate',
            'coral',
            'cornflowerblue',
            #'cornsilk',
            'crimson',
            'cyan',
            'darkblue',
            'darkcyan',
            'darkgoldenrod',
            'darkgreen',
            'darkkhaki',
            'darkmagenta',
            'darkolivegreen',
            'darkorange',
            'darkorchid',
            'darkred',
            'darksalmon',
            'darkseagreen',
            'darkslateblue',
            #'darkslategray',
            'darkturquoise',
            'darkviolet',
            'deeppink',
            'deepskyblue',
            'dodgerblue',
            'firebrick',
            'forestgreen',
            'fuchsia',
            #'gainsboro',
            'gold',
            'goldenrod',
            'green',
            'greenyellow',
            #'honeydew',
            'hotpink',
            'indianred',
            'indigo',
            'khaki',
            #'lavender',
            #'lavenderblush',
            'lawngreen',
            #'lemonchiffon',
            'lightblue',
            #'lightcoral',
            #'lightcyan'
            #'lightgoldenrodyellow',
            'lightgreen',
            'lightpink',
            'lightsalmon',
            'lightseagreen',
            'lightskyblue',
            #'lightslategray',
            #'lightsteelblue',
            #'lightyellow',
            'lime',
            'limegreen',
            'magenta',
            'maroon',
            'mediumaquamarine',
            'mediumblue',
            'mediumorchid',
            'mediumpurple',
            'mediumseagreen',
            'mediumslateblue',
            'mediumspringgreen',
            'mediumturquoise',
            'mediumvioletred',
            'midnightblue',
            #'mintcream',
            #'mistyrose',
            #'moccasin',
            'navy',
            #'oldlace',
            'olive',
            'olivedrab',
            'orange',
            'orangered',
            'orchid',
            #'palegoldenrod',
            'palegreen',
            'paleturquoise',
            'palevioletred',
            #'papayawhip',
            'peachpuff',
            #'peru',
            'pink',
            'plum',
            'powderblue',
            'purple',
            'red',
            'rosybrown',
            'royalblue',
            'saddlebrown',
            'salmon',
            'sandybrown',
            'seagreen',
            #'seashell',
            'sienna',
            'skyblue',
            'slateblue',
            'springgreen',
            'steelblue',
            #'tan',
            'teal',
            #'thistle',
            'tomato',
            'turquoise',
            'violet',
            #'wheat',
            #'yellow',
            #'yellowgreen'
        ]

    def package_folder(self, folder_path, zip_file_name, zip_file_description):
        ''' Simple utility for packaging a folder and saving to shock '''
        if folder_path == self.scratch:
            raise ValueError ("cannot package folder that is not a subfolder of scratch")
        dfu = DataFileUtil(self.callback_url)
        if not os.path.exists(folder_path):
            raise ValueError ("cannot package folder that doesn't exist: "+folder_path)
        output = dfu.file_to_shock({'file_path': folder_path,
                                    'make_handle': 0,
                                    'pack': 'zip'})
        return {'shock_id': output['shock_id'],
                'name': zip_file_name,
                'label': zip_file_description}


    def make_kaiju_summary_plots(self, in_folder_path, bar_plot_out_folder_path, area_plot_out_folder_path):
        ''' Simple utility for packaging a folder and saving to shock '''
        if not os.path.exists(in_folder_path):
            raise ValueError ("kaiju summary folder doesn't exist: "+in_folder_path)

        # read each sample file for each tax level and store abundances
        abundances = dict()

        self._create_bar_plots (abundances)
        self._create_area_plots (abundances)
        


    def build_html_output_for_lineage_wf(self, html_dir, object_name):
        '''
        Based on the output of CheckM lineage_wf, build an HTML report
        '''

        # move plots we need into the html directory
        plot_name = 'bin_qa_plot.png'
        shutil.copy(os.path.join(self.plots_dir, plot_name), os.path.join(html_dir, plot_name))
        self._copy_ref_dist_plots(self.plots_dir, html_dir)

        # write the html report to file
        html = open(os.path.join(html_dir, 'report.html'), 'w')

        # header
        self._write_html_header(html, object_name)
        html.write('<body>\n')

        # include the single main summary figure
        html.write('<img src="' + plot_name + '" width="90%" />\n')
        html.write('<br><br><br>\n')

        # print out the info table
        self.build_summary_table(html, html_dir)

        html.write('</body>\n</html>\n')
        html.close()

        return self.package_folder(html_dir, 'report.html', 'Assembled report from CheckM')


    def build_summary_table(self, html, html_dir):

        stats_file = os.path.join(self.output_dir, 'storage', 'bin_stats_ext.tsv')
        if not os.path.isfile(stats_file):
            log('Warning! no stats file found (looking at: ' + stats_file + ')')
            return

        bin_stats = []
        with open(stats_file) as lf:
            for line in lf:
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                col = line.split('\t')
                bin_id = col[0]
                data = ast.literal_eval(col[1])
                bin_stats.append({'bid': bin_id, 'data': data})


        fields = [{'id': 'marker lineage', 'display': 'Marker Lineage'},
                  {'id': '# genomes', 'display': '# Genomes'},
                  {'id': '# markers', 'display': '# Markers'},
                  {'id': '# marker sets', 'display': '# Marker Sets'},
                  {'id': '0', 'display': '0'},
                  {'id': '1', 'display': '1'},
                  {'id': '2', 'display': '2'},
                  {'id': '3', 'display': '3'},
                  {'id': '4', 'display': '4'},
                  {'id': '5+', 'display': '5+'},
                  {'id': 'Completeness', 'display': 'Completeness', 'round': 3},
                  {'id': 'Contamination', 'display': 'Contamination', 'round': 3}]

        html.write('<table>\n')
        html.write('  <tr>\n')
        html.write('    <th><b>Bin Name</b></th>\n')
        for f in fields:
            html.write('    <th>' + f['display'] + '</th>\n')
        html.write('  </tr>\n')

        for b in bin_stats:
            html.write('  <tr>\n')
            dist_plot_file = os.path.join(html_dir, str(b['bid']) + self.DIST_PLOT_EXT)
            if os.path.isfile(dist_plot_file):
                self._write_dist_html_page(html_dir, b['bid'])
                html.write('    <td><a href="' + b['bid'] + '.html">' + b['bid'] + '</td>\n')
            else:
                html.write('    <td>' + b['bid'] + '</td>\n')
            for f in fields:
                if f['id'] in b['data']:
                    value = str(b['data'][f['id']])
                    if f.get('round'):
                        value = str(round(b['data'][f['id']], f['round']))
                    html.write('    <td>' + value + '</td>\n')
                else:
                    html.write('    <td></td>\n')
            html.write('  </tr>\n')

        html.write('</table>\n')


    def _create_bar_plots (self, abundances):
        color_names = self.no_light_color_names

        import numpy as np
        import matplotlib.pyplot as plt
        import random
        from random import shuffle

        y_label = 'percent'
        title = 'Lineage Proportion'
        sample_labels = ['sample1', 'sample2', 'sample3', 'sample4', 'sample5']
        element_labels = ['OTU_1', 'OTU_2', 'OTU_3', 'OTU_4']

        N = len(sample_labels)
        random.seed(a=len(element_labels))
        r = random.random()
        shuffle(color_names, lambda: r)


        vals = [[20, 35, 20, 35, 27],
                [25, 22, 34, 20, 15],
                [45, 33, 36, 35, 48],
                [10, 10, 10, 10, 10]
            ]
        ind = np.arange(N)    # the x locations for the groups
        bar_width = 0.5      # the width of the bars: can also be len(x) sequence
        label_ind = []
        for ind_i,this_ind in enumerate(ind):
            ind[ind_i] = this_ind+bar_width/2
            label_ind.append(this_ind + bar_width/2)
        np_vals = []
        for vec_i,val_vec in enumerate(vals):
            np_vals.append(np.array(val_vec))
    
        # Build image
        if N < 10:
            img_in_width = 2*N
        elif N < 20:
            img_in_width = N
        else:
            img_in_width = 20
        img_in_height = 5
        fig = plt.figure()
        fig.set_size_inches(img_in_width, img_in_height)
        ax = plt.subplot(111)

        #for ax in fig.axes:
        #    ax.xaxis.set_visible(False)  # remove axis labels and ticks
        #    ax.yaxis.set_visible(False)
        #    for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tick labels
        #        t.set_visible(False)
        #for ax in fig.axes:
        #    ax.spines['top'].set_visible(False) # Get rid of top axis line
        #    ax.spines['bottom'].set_visible(False) #  bottom axis line
        #    ax.spines['left'].set_visible(False) #  Get rid of bottom axis line
        #    ax.spines['right'].set_visible(False) #  Get rid of bottom axis line
    
        last_bottom = None
        p = []
        for vec_i,val_vec in enumerate(np_vals):
            if vec_i == 0:
                this_bottom = 0
                last_bottom = val_vec
            else:
                this_bottom = last_bottom
                last_bottom = this_bottom + val_vec
            p.append (ax.bar (ind, val_vec, bar_width, bottom=this_bottom, color=color_names[vec_i], alpha=0.4, ec='none'))

        plt.ylabel(y_label)
        plt.title(title)
        plt.xticks(label_ind, sample_labels, ha='right', rotation=45)
        plt.yticks(np.arange(0, 101, 10))

        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])    
        key_colors = []
        for each_p in reversed(p):
            key_colors.append(each_p[0])
        plt.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1, 1))

        plt.show()


    def _create_area_plots (self, abundances):
        color_names = self.no_light_color_names
        
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import random
        from random import shuffle
        
        y_label = 'percent'
        title = 'Lineage Proportion'
        sample_labels = ['sample1', 'sample2', 'sample3', 'sample4', 'sample5']
        element_labels = ['OTU_1', 'OTU_2', 'OTU_3', 'OTU_4']

        N = len(sample_labels)
        random.seed(a=len(element_labels))

        r = random.random()
        shuffle(color_names, lambda: r)

        vals = [[20, 35, 20, 35, 27],
                [25, 22, 34, 20, 15],
                [45, 33, 36, 35, 48],
                [10, 10, 10, 10, 10]
            ]
        ind = np.arange(N)    # the x locations for the groups
        label_ind = ind

        np_vals = []
        for vec_i,val_vec in enumerate(vals):
            np_vals.append(np.array(val_vec))
    
        # Build image
        if N < 10:
            img_in_width = 2*N
        elif N < 20:
            img_in_width = N
        else:
            img_in_width = 20
        img_in_height = 5
        fig = plt.figure()
        fig.set_size_inches(img_in_width, img_in_height)
        ax = plt.subplot(111)

        # Let's turn off visibility of all tic labels and boxes here
        #for ax in fig.axes:
        #    ax.xaxis.set_visible(False)  # remove axis labels and tics
        #    ax.yaxis.set_visible(False)
        #    for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tics
        #        t.set_visible(False)
        #    ax.spines['top'].set_visible(False)     # Get rid of top axis line
        #    ax.spines['bottom'].set_visible(False)  # bottom axis line
        #    ax.spines['left'].set_visible(False)    # left axis line
        #    ax.spines['right'].set_visible(False)   # right axis line

        ax.stackplot (ind, np_vals, colors=color_names, alpha=0.4, edgecolor='none')

        plt.ylabel(y_label)
        plt.title(title)
        plt.xticks(label_ind, sample_labels, ha='right', rotation=45)
        plt.yticks(np.arange(0, 101, 10))

        # creating the legend manually
        key_colors = []
        for color_i in reversed(np.arange(N-1)):
            key_colors.append(mpatches.Patch(color=color_names[color_i], alpha=0.4, ec='black'))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])   
        plt.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1, 1))

        #plt.grid()
        plt.show()


    def _write_html_header(self, html, object_name):

        html.write('<html>\n')
        html.write('<head>\n')
        html.write('<title>CheckM Report for ' + object_name + '</title>')

        style = '''
        <style style="text/css">
            a {
                color: #337ab7;
            }

            a:hover {
                color: #23527c;
            }

            table {
                border: 1px solid #bbb;
                border-collapse: collapse;
            }

            th, td {
                text-align: left;
                border: 1px solid #bbb;
                padding: 8px;
            }

            tr:nth-child(odd) {
                background-color: #f9f9f9;
            }

            tr:hover {
                background-color: #f5f5f5;
            }
        </style>\n</head>\n'''

        html.write(style)
        html.write('</head>\n')

    def _copy_file_ignore_errors(self, filename, src_folder, dest_folder):
        src = os.path.join(src_folder, filename)
        dest = os.path.join(dest_folder, filename)
        log('copying ' + src + ' to ' + dest)
        try:
            shutil.copy(src, dest)
        except:
            # TODO: add error message reporting
            log('copy failed')


    def _write_dist_html_page(self, html_dir, bin_id):

        # write the html report to file
        html = open(os.path.join(html_dir, bin_id + '.html'), 'w')

        html.write('<html>\n')
        html.write('<head>\n')
        html.write('<title>CheckM Dist Plots for Bin' + bin_id + '</title>')
        html.write('<style style="text/css">\n a { color: #337ab7; } \n a:hover { color: #23527c; }\n</style>\n')
        html.write('<body>\n')
        html.write('<br><a href="report.html">Back to summary</a><br>\n')
        html.write('<center><h2>Bin: ' + bin_id + '</h2></center>\n')
        html.write('<img src="' + bin_id + self.DIST_PLOT_EXT + '" width="90%" />\n')
        html.write('<br><br><br>\n')
        html.write('</body>\n</html>\n')
        html.close()


    def _copy_ref_dist_plots(self, plots_dir, dest_folder):
        for plotfile in os.listdir(plots_dir):
            plot_file_path = os.path.join(plots_dir, plotfile)
            if os.path.isfile(plot_file_path) and plotfile.endswith(self.DIST_PLOT_EXT):
                try:
                    shutil.copy(os.path.join(plots_dir, plotfile),
                                os.path.join(dest_folder, plotfile))
                except:
                    # TODO: add error message reporting
                    log('copy of ' + plot_file_path + ' to html directory failed')
