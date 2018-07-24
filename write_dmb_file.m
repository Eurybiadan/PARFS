function [  ] = write_dmb_file( dmb_fullpath, dmb_contents )
% [  ] = write_dmb_file( dmb_fullpath, dmb_contents )
% 
% Robert F Cooper, 2018-06-18
%
% This script takes a matlab struct and outputs it to a binary python
% "pickle".


dmb_dict = py.dict(pyargs('frame_strip_ncc_threshold', dmb_contents.frame_strip_ncc_threshold,...
       'n_columns_desinusoided', dmb_contents.n_columns_desinusoided,...       
       'strip_n_frames_with_highest_ncc_value', dmb_contents.strip_n_frames_with_highest_ncc_value,...
       'full_frame_n_frames_with_highest_ncc_value',dmb_contents.full_frame_n_frames_with_highest_ncc_value,...
       'image_sequence_file_name', dmb_contents.image_sequence_file_name,...
       'reference_frame', dmb_contents.reference_frame,...
       'secondary_sequences_file_names', dmb_contents.secondary_sequences_file_names',...
       'secondary_sequences_absolute_paths', dmb_contents.secondary_sequences_absolute_paths',...
       'frame_strip_lines_per_strip', dmb_contents.frame_strip_lines_per_strip,...
       'frame_strip_lines_between_strips_start', dmb_contents.frame_strip_lines_between_strips_start,...
       'n_frames', dmb_contents.n_frames,...
       'save_strip_registered_sequence', dmb_contents.save_strip_registered_sequence,...
       'frame_strip_ncc_n_columns_to_ignore', dmb_contents.frame_strip_ncc_n_columns_to_ignore,...       
       'image_sequence_absolute_path', dmb_contents.image_sequence_absolute_path,...
       'n_columns_raw_sequence', dmb_contents.n_columns_raw_sequence,...
       'fast_scanning_horizontal', dmb_contents.fast_scanning_horizontal,...
       'n_rows_desinusoided', dmb_contents.n_rows_desinusoided,...
       'n_rows_raw_sequence', dmb_contents.n_rows_raw_sequence,...
       'desinusoid_data_filename', dmb_contents.desinusoid_data_filename,...
       'desinusoid_data_absolute_path', dmb_contents.desinusoid_data_absolute_path,...    
       'strip_DCT_terms_retained_percentage', dmb_contents.strip_DCT_terms_retained_percentage,...
       'frame_strip_ncc_n_rows_to_ignore', dmb_contents.frame_strip_ncc_n_rows_to_ignore,...
       'desinusoid_matrix', dmb_contents.desinusoid_matrix,...
       'strip_max_displacement_threshold', dmb_contents.strip_max_displacement_threshold,...
       'full_frame_max_displacement_threshold', dmb_contents.full_frame_max_displacement_threshold,...
       'full_frame_ncc_n_lines_to_ignore', dmb_contents.full_frame_ncc_n_lines_to_ignore,...
       'min_overlap_for_cropping_strip_image', dmb_contents.min_overlap_for_cropping_strip_image,...
       'min_overlap_for_cropping_full_frame_image', dmb_contents.min_overlap_for_cropping_full_frame_image,...
       'strip_registration_required', dmb_contents.strip_registration_required,...
       'save_full_frame_registered_image', dmb_contents.save_full_frame_registered_image,...
       'save_strip_registered_image', dmb_contents.save_strip_registered_image,...
       'frame_strip_calculation_precision', dmb_contents.frame_strip_calculation_precision,...
       'desinusoiding_required', dmb_contents.desinusoiding_required,...
       'clinical_version', dmb_contents.clinical_version,...
       'full_frame_calculation_precision', dmb_contents.full_frame_calculation_precision,...
       'save_full_frame_registered_sequence', dmb_contents.save_full_frame_registered_sequence,...
       'user_defined_suffix', dmb_contents.user_defined_suffix));

thispath = getparent(which('write_dmb_file'));   
   
if count(py.sys.path,thispath) == 0
    insert(py.sys.path,int32(0),thispath);
end

if isunix
    RTLD_NOW = 2;
    RTLD_DEEPBIND = 8;
    flag = bitor(RTLD_NOW, RTLD_DEEPBIND);    
    py.sys.setdlopenflags(int32(flag));
end

py.matlab_pickler.pickle_mat(dmb_fullpath,dmb_dict);   

end

