# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import json

from kb_kaiju.Utils.KaijuUtil import KaijuUtil
#END_HEADER


class kb_kaiju:
    '''
    Module Name:
    kb_kaiju

    Module Description:
    A KBase module: kb_kaiju
This SDK module is developed to wrap the open source package Kaiju
determines taxonomic structure for microbial communities from shotgun metagenomic sequence data

Module also utilizes Krona for visualization of results

References: 
Kaiju Homepage: http://kaiju.binf.ku.dk/
Krona Homepage: https://github.com/marbl/Krona/wiki
Kaiju DBs from: http://kaiju.binf.ku.dk/server
Github repo for Kaiju: https://github.com/bioinformatics-centre/kaiju
Github repo for Krona: https://github.com/marbl/Krona

Kaiju paper: Menzel, P. et al. (2016) Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nat. Commun. 7:11257.
Krona paper: Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.

Kaiju License:
Copyright (c) 2015, 2016, 2017 Peter Menzel and Anders Krogh

Kaiju is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Kaiju is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the file LICENSE for more details.

You should have received a copy of the GNU General Public License along with the source code. If not, see http://www.gnu.org/licenses/.
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.1.0"
    GIT_URL = "https://github.com/kbaseapps/kb_kaiju"
    GIT_COMMIT_HASH = "6d7e531782136056d61638c13f43afc1c9e2989c"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        #END_CONSTRUCTOR
        pass


    def run_kaiju(self, ctx, params):
        """
        Kaiju Method
        :param params: instance of type "KaijuInputParams" (Kaiju App Input
           Params) -> structure: parameter "workspace_name" of type
           "workspace_name" (** The workspace object refs are of form: ** ** 
           objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_refs" of type "data_obj_ref", parameter
           "tax_levels" of list of String, parameter "db_type" of String,
           parameter "filter_percent" of Double, parameter
           "subsample_percent" of Long, parameter "subsample_replicates" of
           Long, parameter "subsample_seed" of Long, parameter "seg_filter"
           of type "bool" (A boolean - 0 for false, 1 for true. @range (0,
           1)), parameter "min_match_length" of Long, parameter
           "greedy_run_mode" of type "bool" (A boolean - 0 for false, 1 for
           true. @range (0, 1)), parameter "greedy_allowed_mismatches" of
           Long, parameter "greedy_min_match_score" of Long, parameter
           "greedy_max_e_value" of Double, parameter "sort_taxa_by" of String
        :returns: instance of type "KaijuOutput" (Kaiju App Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_kaiju
        print('--->\nRunning kb_kaiju.run_kaiju\nparams:')
        print(json.dumps(params, indent=1))

        ku = KaijuUtil(self.config, ctx)
        returnVal = ku.run_kaiju_with_krona(params)
        #END run_kaiju

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_kaiju return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
