
import numpy
import os,pickle


def pickle_mat(pickle_path, this_dict):
	

    this_dict['desinusoid_matrix']=numpy.reshape(this_dict['desinusoid_matrix'], (int(this_dict['n_columns_desinusoided']), int(this_dict['n_columns_raw_sequence']) ) )
    this_dict['desinusoid_matrix']=numpy.transpose(this_dict['desinusoid_matrix'])

    this_dict['secondary_sequences_absolute_paths']=list(this_dict['secondary_sequences_absolute_paths'])
    this_dict['secondary_sequences_file_names']=list(this_dict['secondary_sequences_file_names'])
    if type(this_dict['strip_n_frames_with_highest_ncc_value']) is not int:
        this_dict['strip_n_frames_with_highest_ncc_value']=list(this_dict['strip_n_frames_with_highest_ncc_value'])

    pickle_file = open(pickle_path, 'wb')
    
    pickle.dump(this_dict,pickle_file)

    pickle_file.close()
