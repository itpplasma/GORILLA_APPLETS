function data=load_copy_transport_coef(data,path_RUN,path_data_plots,gorilla_applets)
    if gorilla_applets.GORILLA_APPLETS_NML.idiffcoef_output == 1
        name=strsplit(gorilla_applets.GORILLA_APPLETS_NML.filename_transp_diff_coef,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(gorilla_applets.GORILLA_APPLETS_NML.filename_transp_diff_coef);
        copyfile([path_RUN,'/',gorilla_applets.GORILLA_APPLETS_NML.filename_transp_diff_coef],path_data_plots);
        
    elseif gorilla_applets.GORILLA_APPLETS_NML.idiffcoef_output == 2
        name=strsplit(gorilla_applets.GORILLA_APPLETS_NML.filename_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(gorilla_applets.GORILLA_APPLETS_NML.filename_delta_s_squared);
        copyfile([path_RUN,'/',gorilla_applets.GORILLA_APPLETS_NML.filename_delta_s_squared],path_data_plots);
    
        name=strsplit(gorilla_applets.GORILLA_APPLETS_NML.filename_std_dvt_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(gorilla_applets.GORILLA_APPLETS_NML.filename_std_dvt_delta_s_squared);
        copyfile([path_RUN,'/',gorilla_applets.GORILLA_APPLETS_NML.filename_std_dvt_delta_s_squared],path_data_plots);
        
    else
        name=strsplit(gorilla_applets.GORILLA_APPLETS_NML.filename_transp_diff_coef,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(gorilla_applets.GORILLA_APPLETS_NML.filename_transp_diff_coef);
        copyfile([path_RUN,'/',gorilla_applets.GORILLA_APPLETS_NML.filename_transp_diff_coef],path_data_plots);
        
        name=strsplit(gorilla_applets.GORILLA_APPLETS_NML.filename_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(gorilla_applets.GORILLA_APPLETS_NML.filename_delta_s_squared);
        copyfile([path_RUN,'/',gorilla_applets.GORILLA_APPLETS_NML.filename_delta_s_squared],path_data_plots);
    
        name=strsplit(gorilla_applets.GORILLA_APPLETS_NML.filename_std_dvt_delta_s_squared,'.dat');
        name=strrep(name,'.','_');
        data.(name{1})=load(gorilla_applets.GORILLA_APPLETS_NML.filename_std_dvt_delta_s_squared);
        copyfile([path_RUN,'/',gorilla_applets.GORILLA_APPLETS_NML.filename_std_dvt_delta_s_squared],path_data_plots);
    end
end
