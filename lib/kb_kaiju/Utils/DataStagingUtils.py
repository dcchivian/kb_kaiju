import os
import time
import shutil
import re
import random
#import subprocess
#import glob

from Workspace.WorkspaceClient import Workspace
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from SetAPI.SetAPIServiceClient import SetAPI


class DataStagingUtils(object):

    def __init__(self, config, ctx):
        self.ctx = ctx
        self.scratch = os.path.abspath(config['scratch'])
        self.ws_url = config['workspace-url']
        self.serviceWizardURL = config['srv-wiz-url']
        self.callbackURL = config['SDK_CALLBACK_URL']
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        self.SE_flag = 'SE'
        self.PE_flag = 'PE'

        SERVICE_VER = 'release'

        # readsUtils_Client
        try:
            self.readsUtils_Client = ReadsUtils(self.callbackURL, token=self.ctx['token'], service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate readsUtils_Client with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))

        # setAPI_Client
        try:
            #setAPI_Client = SetAPI (url=self.callbackURL, token=self.ctx['token'])  # for SDK local.  local doesn't work for SetAPI
            self.setAPI_Client = SetAPI (url=self.serviceWizardURL, token=self.ctx['token'])  # for dynamic service
        except Exception as e:
            raise ValueError('Unable to instantiate setAPI_Client with serviceWizardURL: '+ self.serviceWizardURL +' ERROR: ' + str(e))


    def expand_input(self, input_refs):
        '''
        Expand input based on an input data reference for Kaiju

        input_refs can be a list of references to a PairedEndLibrary, a SingleEndLibrary, or a ReadsSet
        '''
        # config
        #SERVICE_VER = 'dev'
        SERVICE_VER = 'release'

        # expand any sets and build a non-redundant list of reads input objs
        ws = Workspace(self.ws_url)
        expanded_input = []
        input_ref_seen = dict()
        SE_types = ['KBaseFile.SingleEndLibrary', 'KBaseAssembly.SingleEndLibrary']
        PE_types = ['KBaseFile.PairedEndLibrary', 'KBaseAssembly.PairedEndLibrary']

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        for input_ref in input_refs:
            input_info = ws.get_object_info3({'objects': [{'ref': input_ref}]})['infos'][0]
            obj_name = input_info[NAME_I]
            type_name = input_info[TYPE_I].split('-')[0]

            # ReadsSet
            if type_name in ['KBaseSets.ReadsSet']:
                try:
                    input_readsSet_obj = self.setAPI_Client.get_reads_set_v1 ({'ref':input_ref,'include_item_info':1})

                except Exception as e:
                    raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(input_ref)+")\n" + str(e))

                for readsLibrary_obj in input_readsSet_obj['data']['items']:
                    this_reads_ref = readsLibrary_obj['ref']
                    if this_reads_ref in input_ref_seen:
                        continue
                    input_ref_seen[this_reads_ref] = True

                    this_reads_name = readsLibrary_obj['info'][NAME_I]
                    reads_item_type = readsLibrary_obj['info'][TYPE_I]
                    reads_item_type = re.sub ('-[0-9]+\.[0-9]+$', "", reads_item_type)  # remove trailing version
                    if reads_item_type in PE_types:
                        this_reads_type = self.PE_flag
                    elif reads_item_type in SE_types:
                        this_reads_type = self.SE_flag
                    else:
                        raise ValueError ("Can't handle read item type '"+reads_item_type+"' obj_name: '"+this_reads_name+" in Set: '"+str(input_ref)+"'")
                    expanded_input.append({'ref':  this_reads_ref,
                                           'name': this_reads_name,
                                           'type': this_reads_type
                                       })
            # SingleEnd Library
            elif type_name in SE_types:
                this_reads_ref = input_ref
                if this_reads_ref in input_ref_seen:
                    continue
                input_ref_seen[this_reads_ref] = True
                this_reads_name = obj_name
                this_reads_type = self.SE_flag
                expanded_input.append({'ref':  this_reads_ref,
                                       'name': this_reads_name,
                                       'type': this_reads_type
                                   })
            # PairedEnd Library
            elif type_name in PE_types:
                this_reads_ref = input_ref
                if this_reads_ref in input_ref_seen:
                    continue
                input_ref_seen[this_reads_ref] = True
                this_reads_name = obj_name
                this_reads_type = self.PE_flag
                expanded_input.append({'ref':  this_reads_ref,
                                       'name': this_reads_name,
                                       'type': this_reads_type
                                   })
            else:
                raise ValueError ("Illegal type in input_refs: "+str(obj_name)+" ("+str(input_ref)+") is of type: '"+str(type_name)+"'")


        return expanded_input


    def stage_input(self,
                    input_item=None,
                    subsample_percent=10,
                    subsample_replicates=1,
                    subsample_seed=1,
                    fasta_file_extension='fastq'):
        '''
        Stage input based on an input data reference for Kaiju

        input_refs can be a list of references to a PairedEndLibrary, a SingleEndLibrary, or a ReadsSet

        This method creates a directory in the scratch area with the set of Fasta/Fastq files, names
        will have the fasta_file_extension parameter tacked on.

            ex:

            staged_input = stage_input({'ref':<ref>,'name':<name>,'type':<type>}, subsample_percent, subsample_replicates, subsample_seed, 'fastq')

            staged_input
            {"input_dir": '...'}
        '''
        # init
        staged_input = dict()
        replicate_input = []

        # config
        #SERVICE_VER = 'dev'
        SERVICE_VER = 'release'

        # generate a folder in scratch to hold the input
        suffix = str(int(time.time() * 1000))
        input_dir = os.path.join(self.scratch, 'input_reads_' + suffix)
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

        #
        # Download reads
        #

        # Paired End Lib
        if input_item['type'] == self.PE_flag:
            try:
                readsLibrary = self.readsUtils_Client.download_reads ({'read_libraries': [input_item['ref']],
                                                                       'interleaved': 'false'})
            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(input_item['ref']) +")\n" + str(e))

            input_fwd_file_path = readsLibrary['files'][input_item['ref']]['files']['fwd']
            input_rev_file_path = readsLibrary['files'][input_item['ref']]['files']['rev']
            fwd_filename = os.path.join(input_dir, input_item['name'] + '.fwd.' + fasta_file_extension)
            rev_filename = os.path.join(input_dir, input_item['name'] + '.rev.' + fasta_file_extension)
            if input_fwd_file_path != fwd_filename:
                shutil.move(input_fwd_file_path, fwd_filename)
            if input_rev_file_path != rev_filename:
                shutil.move(input_rev_file_path, rev_filename)
            input_item['fwd_file'] = fwd_filename
            input_item['rev_file'] = rev_filename

            if not os.path.isfile(fwd_filename):
                raise ValueError('Error generating reads file '+fwd_filename)
            if not os.path.isfile(rev_filename):
                raise ValueError('Error generating reads file '+rev_filename)
            # make sure fasta file isn't empty
            min_fasta_len = 1
            if not self._fasta_seq_len_at_least(fwd_filename, min_fasta_len):
                raise ValueError('Reads Library is empty in filename: '+str(fwd_filename))
            if not self._fasta_seq_len_at_least(rev_filename, min_fasta_len):
                raise ValueError('Reads Library is empty in filename: '+str(rev_filename))

        # Single End Lib
        elif input_item['type'] == self.SE_flag:
            try:
                readsLibrary = self.readsUtils_Client.download_reads ({'read_libraries': [input_item['ref']]})
            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(input_item['ref']) +")\n" + str(e))

            input_fwd_file_path = readsLibrary['files'][input_item['ref']]['files']['fwd']
            fwd_filename = os.path.join(input_dir, input_item['name'] + '.fwd.' + fasta_file_extension)
            if input_fwd_file_path != fwd_filename:
                shutil.move(input_fwd_file_path, fwd_filename)
            input_item['fwd_file'] = fwd_filename

            if not os.path.isfile(fwd_filename):
                raise ValueError('Error generating reads file '+fwd_filename)
            # make sure fasta file isn't empty
            min_fasta_len = 1
            if not self._fasta_seq_len_at_least(fwd_filename, min_fasta_len):
                raise ValueError('Reads Library is empty in filename: '+str(fwd_filename))

        else:
            raise ValueError ("No type set for input library "+str(input_item['name'])+" ("+str(input_item['ref'])+")")


        #
        # Subsample
        #

        if subsample_percent == 100:
            replicate_input = [input_item]
        else:
            replicate_input = self._randomly_subsample_reads(input_item,
                                                             subsample_percent    = subsample_percent,
                                                             subsample_replicates = subsample_replicates,
                                                             subsample_seed       = subsample_seed)
            # free up disk
            os.remove(input_item['fwd_file'])
            if input_item['type'] == self.PE_flag:
                os.remove(input_item['rev_file'])


        # return input file info
        #staged_input['input_dir'] = input_dir
        #staged_input['folder_suffix'] = suffix
        staged_input['replicate_input'] = replicate_input
        return staged_input


    def _randomly_subsample_reads(self,
                                  input_item=None,
                                  subsample_percent=100,
                                  subsample_replicates=1,
                                  subsample_seed=1):

        replicate_files = []
        split_num = subsample_replicates

        # for now can only do percentage instead of raw cnt of reads per subsample
        use_reads_num  = False
        use_reads_perc = True
        reads_num = 0  # not used.  subsample_percent used instead

        # init randomizer
        random.seed(subsample_seed)


        # Paired End
        #
        if input_item['type'] == self.PE_flag:
            print ("SUBSAMPLING PE library "+input_item['name'])  # DEBUG

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_item['fwd_file'])
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            input_rev_path = re.sub ("\.fastq$", "", input_item['rev_file'])
            input_rev_path = re.sub ("\.FASTQ$", "", input_rev_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"
            output_rev_paired_file_path_base   = input_rev_path+"_rev_paired"

            # set up for file io
            total_paired_reads = 0
            total_unpaired_fwd_reads = 0
            total_unpaired_rev_reads = 0
            total_paired_reads_by_set = []
            fwd_ids = dict()
            paired_ids = dict()
            paired_ids_list = []
            paired_lib_i = dict()
            paired_buf_size = 100000
            recs_beep_n = 1000000


            # read fwd file to get fwd ids
#            rec_cnt = 0  # DEBUG
            print ("GETTING IDS")  # DEBUG
            with open (input_item['fwd_file'], 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        # added below line to manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        fwd_ids[read_id] = True

                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            print ("read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1


            # read reverse to determine paired
            print ("DETERMINING PAIRED IDS")  # DEBUG
            with open (input_item['rev_file'], 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        # added below line to manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        if fwd_ids.get(read_id, False):
                            paired_ids[read_id] = True
                            paired_ids_list.append(read_id)

                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            print ("read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1
            total_paired_reads = len(paired_ids_list)
            print ("TOTAL PAIRED READS CNT: "+str(total_paired_reads))  # DEBUG


            # Determine sublibrary sizes
            if use_reads_num:
                reads_per_lib = reads_num
                if reads_per_lib > total_paired_reads // split_num:
                    raise ValueError ("must specify reads_num <= total_paired_reads_cnt / split_num.  You have reads_num:"+str(reads_num)+" > total_paired_reads_cnt:"+str(total_paired_reads)+" / split_num:"+str(split_num)+".  Instead try reads_num <= "+str(total_paired_reads // split_num))
            elif use_reads_perc:
                reads_per_lib = int ((subsample_percent/100.0) * total_paired_reads)
                if reads_per_lib > total_paired_reads // split_num:
                    raise ValueError ("must specify reads_perc <= 1 / split_num.  You have reads_perc:"+str(subsample_percent)+" > 1 / split_num:"+str(split_num)+".  Instead try reads_perc <= "+ str(int(100 * 1/split_num)))
            else:
                raise ValueError ("error in logic reads_num vs. reads_perc logic")


            # Determine random membership in each sublibrary
            print ("GETTING RANDOM SUBSAMPLES")  # DEBUG
            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * split_num)):
                lib_i = i % split_num
                paired_lib_i[read_id] = lib_i


            # split fwd paired
            print ("WRITING FWD SPLIT PAIRED")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(split_num):
                paired_output_reads_file_handles.append(open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))
                total_paired_reads_by_set.append(0)

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_item['fwd_file'], 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_lib_i[last_read_id]
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                                total_paired_reads_by_set[lib_i] += 1
                                if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                    print ("\t"+str(paired_cnt)+" recs processed")
                            else:
                                #unpaired_fwd_buf.extend(rec_buf)
                                pass
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        # added below line to manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
                            found = paired_lib_i[read_id]
                            capture_type_paired = True
                        except:
                            total_unpaired_fwd_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if capture_type_paired:
                        lib_i = paired_lib_i[last_read_id]
                        paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                        paired_cnt += 1
                        if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            print ("\t"+str(paired_cnt)+" recs processed")
                    else:
                        #unpaired_fwd_buf.extend(rec_buf)
                        pass
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            print ("\t"+str(paired_cnt)+" FWD recs processed")


            # split rev paired
            print ("WRITING REV SPLIT PAIRED")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(split_num):
                paired_output_reads_file_handles.append(open (output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_item['rev_file'], 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_lib_i[last_read_id]
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                                if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                    print ("\t"+str(paired_cnt)+" recs processed")
                            else:
                                #unpaired_fwd_buf.extend(rec_buf)
                                pass
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        # added below line to manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
                            found = paired_lib_i[read_id]
                            capture_type_paired = True
                        except:
                            total_unpaired_rev_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if capture_type_paired:
                        lib_i = paired_lib_i[last_read_id]
                        paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                        paired_cnt += 1
                        if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            print ("\t"+str(paired_cnt)+" recs processed")
                    else:
                        #unpaired_fwd_buf.extend(rec_buf)
                        pass
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            print ("\t"+str(paired_cnt)+" REV recs processed")


            # summary
            report = 'SUMMARY FOR SUBSAMPLE OF READ LIBRARY: '+input_item['name']+"\n"
            report += "TOTAL PAIRED READS: "+str(total_paired_reads)+"\n"
            report += "TOTAL UNPAIRED FWD READS (discarded): "+str(total_unpaired_fwd_reads)+"\n"
            report += "TOTAL UNPAIRED REV READS (discarded): "+str(total_unpaired_rev_reads)+"\n"
            report += "\n"
            for lib_i in range(split_num):
                report += "PAIRED READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"
            print (report)


            # make replicate objects to return
#        for replicate_i,replicate_item in enumerate(replicate_files):
#            replicate_input.append({'fwd_file': replicate_item['fwd_file'],
#                                    'type': input_item['type'],
#                                    'name': input_item['name']+"-"+str(replicate_i)
#                                })
#            if input_item['type'] == self.PE_flag:
#                replicate_input[replicate_i]['rev_file'] = replicate_item['rev_file']

            print ("MAKING REPLICATE OBJECT")  # DEBUG
            paired_obj_refs = []
            for lib_i in range(split_num):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                output_rev_paired_file_path = output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                     or os.path.getsize (output_fwd_paired_file_path) == 0 \
                   or not os.path.isfile (output_rev_paired_file_path) \
                     or os.path.getsize (output_rev_paired_file_path) == 0:

                    raise ValueError ("failed to create paired output")
                else:
                    zero_pad = '0'*(len(str(split_num))-len(str(lib_i+1)))
                    replicate_files.append({'fwd_file': output_fwd_paired_file_path,
                                            'rev_file': output_rev_paired_file_path,
                                            'ref':  input_item['ref'],  # note: this is for the src, not the subsample which is not saved
                                            'type': input_item['type'],
                                            'name': input_item['name']+'-'+zero_pad+str(lib_i+1)
                                        })

        # SingleEndLibrary
        #
        elif input_item['type'] == self.SE_flag:
            print ("SUBSAMPLING SE library "+input_item['name'])

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_item['fwd_file'])
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"

            # get "paired" ids
            print ("DETERMINING IDS")  # DEBUG
            paired_ids = dict()
            paired_ids_list = []
            paired_lib_i = dict()
            paired_buf_size = 100000
            recs_beep_n = 100000

            with open (input_item['fwd_file'], 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        # added below line to manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        if read_id in paired_ids:
                            raise ValueError ("repeat read_id: "+read_id)
                        paired_ids[read_id] = True
                        paired_ids_list.append(read_id)
                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            print ("read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1
            total_paired_reads = len(paired_ids_list)
            print ("TOTAL READS CNT: "+str(total_paired_reads))  # DEBUG


            # Determine sublibrary sizes
            if use_reads_num:
                reads_per_lib = reads_num
                if reads_per_lib > total_paired_reads // split_num:
                    raise ValueError ("must specify reads_num <= total_paired_reads_cnt / split_num.  You have reads_num:"+str(reads_num)+" > total_paired_reads_cnt:"+str(total_paired_reads)+" / split_num:"+str(split_num)+".  Instead try reads_num <= "+str(total_paired_reads // split_num))
            elif use_reads_perc:
                reads_per_lib = int ((subsample_percent/100.0) * total_paired_reads)
                if reads_per_lib > total_paired_reads // split_num:
                    raise ValueError ("must specify reads_perc <= 1 / split_num.  You have reads_perc:"+str(subsample_percent)+" > 1 / split_num:"+str(split_num)+".  Instead try reads_perc <= "+ str(int(100 * 1/split_num)))
            else:
                raise ValueError ("error in logic reads_num vs. reads_perc logic")


            # Determine random membership in each sublibrary
            print ("GETTING RANDOM SUBSAMPLES")  # DEBUG
            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * split_num)):
                lib_i = i % split_num
                paired_lib_i[read_id] = lib_i


            # set up for file io
            total_paired_reads = 0
            total_paired_reads_by_set = []
            paired_buf_size = 1000000


            # split reads
            print ("WRITING SPLIT SINGLE END READS")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(split_num):
                paired_output_reads_file_handles.append(open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))
                total_paired_reads_by_set.append(0)

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            recs_beep_n = 1000000
            with open (input_item['fwd_file'], 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        total_paired_reads += 1
                        if last_read_id != None:
                            try:
                                lib_i = paired_lib_i[last_read_id]
                                total_paired_reads_by_set[lib_i] += 1
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                            except:
                                pass
                            if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                print ("\t"+str(paired_cnt)+" recs processed")
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        # added below line to manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if last_read_id != None:
                        try:
                            lib_i = paired_lib_i[last_read_id]
                            total_paired_reads_by_set[lib_i] += 1
                            paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                            paired_cnt += 1
                        except:
                            pass
                    if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                        print ("\t"+str(paired_cnt)+" recs processed")
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            # summary
            report = 'SUMMARY FOR SUBSAMPLE OF READ LIBRARY: '+input_item['name']+"\n"
            report += "TOTAL READS: "+str(total_paired_reads)+"\n"
            for lib_i in range(split_num):
                report += "SINGLE END READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"
            print (report)


            # make replicate objects to return
            print ("MAKING REPLICATE OBJECTS")  # DEBUG
            paired_obj_refs = []
            for lib_i in range(split_num):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                     or os.path.getsize (output_fwd_paired_file_path) == 0:

                    raise ValueError ("failed to create paired output")
                else:
                    zero_pad = '0'*(len(str(split_num))-len(str(lib_i+1)))
                    replicate_files.append({'fwd_file': output_fwd_paired_file_path,
                                            'ref':  input_item['ref'],  # note: this is for the src, not the subsample which is not saved
                                            'type': input_item['type'],
                                            'name': input_item['name']+'-'+zero_pad+str(lib_i+1)
                                        })


        else:
            raise ValueError ("unknown ReadLibrary type:"+str(input_item['type'])+" for readslibrary: "+input_item['name'])


        return replicate_files


    def _fasta_seq_len_at_least(self, fasta_path, min_fasta_len=1):
        '''
        counts the number of non-header, non-whitespace characters in a FASTA file
        '''
        seq_len = 0
        with open (fasta_path, 'r', 0) as fasta_handle:
            for line in fasta_handle:
                line = line.strip()
                if line.startswith('>'):
                    continue
                line = line.replace(' ','')
                seq_len += len(line)
                if seq_len >= min_fasta_len:
                    return True
        return False
