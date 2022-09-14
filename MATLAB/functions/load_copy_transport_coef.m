function data=load_copy_transport_coef(data,path_RUN,path_data_plots,mono_energetic_transp_coef)
    if mono_energetic_transp_coef.TRANSPCOEFNML.idiffcoef_output == 1
        name=strsplit(mono_energetic_transp_coef.TRANSPCOEFNML.filename_transp_diff_coef,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(mono_energetic_transp_coef.TRANSPCOEFNML.filename_transp_diff_coef);
        copyfile([path_RUN,'/',mono_energetic_transp_coef.TRANSPCOEFNML.filename_transp_diff_coef],path_data_plots);
        
    elseif mono_energetic_transp_coef.TRANSPCOEFNML.idiffcoef_output == 2
        name=strsplit(mono_energetic_transp_coef.TRANSPCOEFNML.filename_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(mono_energetic_transp_coef.TRANSPCOEFNML.filename_delta_s_squared);
        copyfile([path_RUN,'/',mono_energetic_transp_coef.TRANSPCOEFNML.filename_delta_s_squared],path_data_plots);
    
        name=strsplit(mono_energetic_transp_coef.TRANSPCOEFNML.filename_std_dvt_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(mono_energetic_transp_coef.TRANSPCOEFNML.filename_std_dvt_delta_s_squared);
        copyfile([path_RUN,'/',mono_energetic_transp_coef.TRANSPCOEFNML.filename_std_dvt_delta_s_squared],path_data_plots);
        
    else
        name=strsplit(mono_energetic_transp_coef.TRANSPCOEFNML.filename_transp_diff_coef,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(mono_energetic_transp_coef.TRANSPCOEFNML.filename_transp_diff_coef);
        copyfile([path_RUN,'/',mono_energetic_transp_coef.TRANSPCOEFNML.filename_transp_diff_coef],path_data_plots);
        
        name=strsplit(mono_energetic_transp_coef.TRANSPCOEFNML.filename_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(mono_energetic_transp_coef.TRANSPCOEFNML.filename_delta_s_squared);
        copyfile([path_RUN,'/',mono_energetic_transp_coef.TRANSPCOEFNML.filename_delta_s_squared],path_data_plots);
    
        name=strsplit(mono_energetic_transp_coef.TRANSPCOEFNML.filename_std_dvt_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(mono_energetic_transp_coef.TRANSPCOEFNML.filename_std_dvt_delta_s_squared);
        copyfile([path_RUN,'/',mono_energetic_transp_coef.TRANSPCOEFNML.filename_std_dvt_delta_s_squared],path_data_plots);
    end
end
