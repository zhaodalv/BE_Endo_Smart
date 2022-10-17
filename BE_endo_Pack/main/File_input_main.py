#!/bin/env python
import sys
import pandas as pd
import numpy as np
sys.path.append("..")
import BE_endo_Pack.Preprocess.Preprocessing as Prepro
import BE_endo_Pack.Effiency.Eff as Effiency
import BE_endo_Pack.Proportion.proportion as Proportion
def pip_file(file_in,Editor,sgRNA_bed_path,workdir,out_file):
   model_path=""
   pre_path=""
   if Editor=='CBE':
      model_path=workdir+"/BE_endo_smart/model_data/CBE_all_420.h5"
      pre_path=workdir+'/BE_endo_smart/model_data/CBE_proportion_pre.npy'
   else:
      model_path=workdir+"/BE_endo_smart/model_data/ABE_H3K27ac_220.h5" 
      pre_path=workdir+'/BE_endo_smart/model_data/ABE_proportion_pre.npy'
   Input_process_obj=Prepro.input_process(file_in,Editor) ##create preprocess obj 
   seq_40s,targets,labels,positions=Input_process_obj.process() ## obtain seq_40s ,targets,labels,positions
   position_df=pd.DataFrame(positions) ##sgRNA position write to local for intersection
   out_file_path=Input_process_obj.write_out(position_df,sgRNA_bed_path) ###write out
   Model_prepare_obj=Prepro.Endo_Prepare(Editor,out_file_path) ##create Model_prepare_obj
   Model_prepare_obj.intersecting(workdir+"/BE_endo_smart/main/BD_intersect.sh",workdir)
   input_df=Model_prepare_obj.Prepare_inputs()
   Model_prepare_obj.read_ins(workdir+"/BE_endo_smart/main/Intersection_temp/result.bed")
   merge_result=Model_prepare_obj.factor_process(input_df)
   ###prediction
   Seq_obj=Effiency.Sequence_prepare(seq_40s)
   X_array=Seq_obj.one_hot_encoding()
   eff_prediction=Effiency.pip_endo(model_path,seq_40s,merge_result.iloc[:,3:].values,labels)
   eff_table=pd.DataFrame(eff_prediction)
   eff_table.columns=['base1', 'base2', 'base3', 'base4', 'base5', 'base6', 'base7', 'base8',
   'base9', 'base10', 'base11', 'base12', 'base13', 'base14', 'base15',
   'base16', 'base17', 'base18', 'base19', 'base20']
   eff_table['seq']=seq_40s
   eff_table['chr']=input_df['sgchrom'].values
   eff_table['start']=input_df['sgstart'].values
   eff_table['end']=input_df['sgend'].values
   eff_table[['chr','start','end','seq','base1', 'base2', 'base3', 'base4', 'base5', 'base6', 'base7', 'base8',
   'base9', 'base10', 'base11', 'base12', 'base13', 'base14', 'base15',
   'base16', 'base17', 'base18', 'base19', 'base20']].to_csv(out_file+".{}.eff".format(Editor),sep="\t",index=None)
   proportion_pre=np.load(pre_path)
   Proportion_prediction_obj=Proportion.Proportion_prediction(Editor,proportion_pre)
   Proportion_df=Proportion_prediction_obj.prediction(seq_40s,targets,labels,eff_prediction)
   Proportion_df.to_csv(out_file+".{}.proportion".format(Editor),sep="\t",index=None)
