
% This script is an moderately modified implementation of the algorithm outlined by Salmon et
% al 2017: "An Automated Reference Frame Selection (ARFS) Algorithm for 
% Cone Imaging with Adaptive Optics Scanning Light Ophthalmoscopy".
%
% At UPenn, we have decided to call it "PARFS: Pretty Accurate Reference
% Frame Selection"
%
% Copyright (C) 2018  Robert F Cooper, Created 06-30-2017
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program. If not, see <http://www.gnu.org/licenses/>.

clear;
close all force;

% Check this version of PARFS against git.
fid = fopen(fullfile(getparent(which(mfilename)),'.VERSION'),'r');
if fid ~= -1
    thisver = fscanf(fid,'%s');
    fclose(fid);
    
    git_version_check( 'Eurybiadan', 'PARFS', thisver )
        
end

params = PARF_Params;

NUM_REF_OUTPUT = params.NUM_REF_OUTPUT; % 5;
MODALITIES = params.MODALITIES; % {'confocal','split_det','avg'}; % The modalities to search for.
MODALITY_WEIGHTS = 1-params.MODALITY_WEIGHTS; % [1/3 1/3 1/3]; % The weights applied to each modality. Adjust these values if you want the first modality (say, confocal) to carry more weight in the reference frame choice.
STRIP_SIZE = params.STRIP_SIZE; % 40; % The size of the strip at which we'll analyze the distortion
BAD_STRIP_THRESHOLD = params.BAD_STRIP_THRESHOLD; % 0; % Having more bad strips than this will result in a frame's removal from consideration.
MIN_NUM_FRAMES_PER_GROUP = params.MIN_NUM_FRAMES_PER_GROUP; % 5; % A group must have more than this number of frames otherwise it will be dropped from consideration

LPS = int32(params.LPS); % int32(12);
LBSS = int32(params.LBSS); % int32(6);
OVERLAP = int32(params.OVERLAP); % int32(5);
NUM_FRAMES = int32(params.NUM_FRAMES); % int32(50);
THRESHOLD = params.THRESHOLD; % 0.7;
OUTPUT_TIFS = params.OUTPUT_TIFS; % true;
OUTPUT_FF_TIFS = params.OUTPUT_FF_TIFS; % false;
OUTPUT_AVIS = params.OUTPUT_AVIS; % false;


fNames=[];
mov_path=[];
addanother=1;

sel_path=pwd;
while addanother ~= 0

    sel_path = uigetdir(sel_path,'Select the folder containing the movies you wish to examine:');

    if sel_path == 0
       break; 
    end
    
    names = read_folder_contents(sel_path,'avi');
    
    for f=1:length(names)
        names{f} = fullfile(sel_path, names{f});        
    end
    
    fNames = [fNames; names];
    mov_path = [mov_path; {sel_path}];
    
    
    choice = questdlg('Would you like to select another folder?','','Yes!','No, thank you.','Yes!');
    
    addanother = strcmp(choice,'Yes!');
end
clear names;

found = zeros(length(MODALITIES),1);
dataSummary = cell(1,1);
dataSummary{1} = 'Input Data Summary:';

stack_fname = cell(length(fNames), length(MODALITIES));
for f=1:length(fNames)
    if ~isempty(fNames{f})
        
        for m = 1:length(MODALITIES) % Find where the modality string is- everything after that, consider to be unique.
            modeind = strfind(fNames{f}, MODALITIES{m});
            if ~isempty(modeind)
               modeind = modeind(end)+length(MODALITIES{m});
               break; 
            end
        end
        
        searchind = fNames{f}(modeind:end);
        samevideos = sort(fNames(~cellfun(@isempty, strfind(fNames, searchind)))); 

        for m = 1:length(MODALITIES)

            matchmode= ~cellfun(@isempty, strfind(samevideos, MODALITIES{m}));
            if any(matchmode)
                stack_fname{f,m} = samevideos{matchmode};
                % Remove it from consideration if we found a match.
                fNames{~cellfun(@isempty, strfind(fNames, samevideos{matchmode}))} = '';
            end

        end
    end
end

for i=1:size(stack_fname,1)
    keep_row(i) = any(~cellfun(@isempty,stack_fname(i,:)));    
end
stack_fname = stack_fname(keep_row,:);


[lut_fname, lut_path]=uigetfile({ '*.csv', 'LUT Files (*.csv)'; '*.xlsx', 'Legacy LUT Format (*.xlsx)' }, ...
                                'Select LUT file:', mov_path{1});
if strcmp(lut_fname(end-3:end), '.csv')
    lut = readtable(fullfile(lut_path,lut_fname));
    headers = lut.Properties.VariableNames;
    version = headers{1};

    switch(version)            
        case 'v0_1'
            lut = table2cell(lut);
            lut(:,1) = cellfun(@(x) ['_' num2str( x,'%04.0f') '_'], lut(:,1),'UniformOutput',false);
            scales = cell2mat(lut(:,4));
    end
else
    version = 'legacy';
    [~,~,lut]=xlsread(fullfile(lut_path,lut_fname));
    scales = cell2mat(lut(:,4));
end

if all(ischar(scales))
    scales = str2num(scales);
end

if ~exist('contains','builtin')
    contains = @(s,p) ~isempty(strfind(s,p));
end

keep_row = false(size(lut,1),1);
for l=1:size(lut,1)
    for f=1 : size(stack_fname,1)
        if ~all(isnan(lut{l,1})) && contains(stack_fname{f,1},lut{l,1})
            keep_row(l)=true;
            continue;
        end
    end
end
lut=lut(keep_row,:);


unique_scales = unique(scales);
unique_scales(isnan(unique_scales)) = [];

dmb_file_to_load=cell(size(stack_fname,1),1);
dmb_path_to_load=cell(size(stack_fname,1),1);

missingentries = size(stack_fname,1)-size(lut,1);
if missingentries>0
    warning('There are missing entries in your LUT file! Files without LUT will be ignored!');
end

for p=1:size(unique_scales,1)
    switch(version)            
        case 'v0_1'
            [dmb_fname, dmb_path]=uigetfile(fullfile(lut_path,'*.mat'),['Select the ***' num2str(unique_scales(p)) '*** degree FOV DESINUSOIDING file!' ]);
        case 'legacy'
            [dmb_fname, dmb_path]=uigetfile(fullfile(lut_path,'*.mat'),['Select the ***' num2str(unique_scales(p)) '*** pixels per fringe DESINUSOIDING file!' ]);
    end
    
    matching_fringes = padarray(scales==unique_scales(p),[missingentries 0],0,'post');
    dmb_file_to_load( matching_fringes )={dmb_fname};
    dmb_path_to_load( matching_fringes )={dmb_path};
end

load(fullfile(dmb_path_to_load{1}, dmb_file_to_load{1}),'horizontal_fringes_n_rows','vertical_fringes_desinusoid_matrix');

desinusoid_matrix = vertical_fringes_desinusoid_matrix;



default_dmb_contents = struct('frame_strip_ncc_threshold', THRESHOLD,...
       'n_columns_desinusoided', int32(size(vertical_fringes_desinusoid_matrix,1)),...
       'n_columns_raw_sequence',int32(size(vertical_fringes_desinusoid_matrix,2)),...
       'strip_n_frames_with_highest_ncc_value', NUM_FRAMES,...
       'full_frame_n_frames_with_highest_ncc_value',NUM_FRAMES,...
       'image_sequence_file_name', stack_fname{f,1},...
       'reference_frame', 0,...
       'secondary_sequences_file_names', [],...
       'secondary_sequences_absolute_paths', [],...
       'frame_strip_lines_per_strip', LPS,...
       'frame_strip_lines_between_strips_start', LBSS,...
       'n_frames', 0,...
       'save_strip_registered_sequence', OUTPUT_AVIS,...
       'frame_strip_ncc_n_columns_to_ignore', int32(150),...       
       'image_sequence_absolute_path', mov_path{1},...       
       'fast_scanning_horizontal', true,...
       'n_rows_desinusoided', int32(horizontal_fringes_n_rows),...
       'n_rows_raw_sequence', int32(horizontal_fringes_n_rows),...
       'desinusoid_data_filename', dmb_file_to_load{1},...
       'desinusoid_data_absolute_path', dmb_path_to_load{1},...    
       'strip_DCT_terms_retained_percentage', int32(50),...
       'frame_strip_ncc_n_rows_to_ignore', int32(3),...
       'desinusoid_matrix', desinusoid_matrix(:),...
       'strip_max_displacement_threshold', int32(200),...
       'full_frame_max_displacement_threshold', int32(200),...
       'full_frame_ncc_n_lines_to_ignore', int32(150),...
       'min_overlap_for_cropping_strip_image',OVERLAP,...
       'min_overlap_for_cropping_full_frame_image',OVERLAP,...
       'strip_registration_required', true,...
       'save_full_frame_registered_image', OUTPUT_FF_TIFS,...
       'save_strip_registered_image', OUTPUT_TIFS,...
       'frame_strip_calculation_precision', 'single',...
       'desinusoiding_required', true,...
       'clinical_version', false,...
       'full_frame_calculation_precision', 'single',...
       'save_full_frame_registered_sequence', false,...
       'user_defined_suffix', ['_ref_0_lps_' num2str(LPS) '_lbss_' num2str(LBSS) ]);
   
% [stack_fname, mov_path] = uigetfile(fullfile(dmb_path,'*.avi'),'Select movies from single timepoint:', 'MultiSelect', 'on');

% stack_fname = {stack_fname};
% mov_path = mov_path;
h=waitbar(0,'Finding reference frames...');

bestrefs=[];

delete(fullfile(mov_path{1},'Reference_Frames.csv'));

%% Analyze the list.
refs = cell(size(stack_fname));
num_frames = zeros(size(stack_fname));

for f=1 : size(stack_fname,1)
    if ~isempty(dmb_path_to_load{f})
        load(fullfile(dmb_path_to_load{f}, dmb_file_to_load{f}),'horizontal_fringes_n_rows','vertical_fringes_desinusoid_matrix');

        default_dmb_contents.desinusoid_matrix = vertical_fringes_desinusoid_matrix';
        default_dmb_contents.desinusoid_matrix = default_dmb_contents.desinusoid_matrix(:)';
        default_dmb_contents.n_rows_desinusoided = int32(horizontal_fringes_n_rows);
        default_dmb_contents.n_rows_raw_sequence = int32(horizontal_fringes_n_rows);
        default_dmb_contents.n_columns_desinusoided= int32(size(vertical_fringes_desinusoid_matrix,1));
        default_dmb_contents.n_columns_raw_sequence= int32(size(vertical_fringes_desinusoid_matrix,2));

        for m=1:size(stack_fname,2)
            if ~isempty(stack_fname{f,m})
                waitbar(f/size(stack_fname,1),h,['Processing video #' num2str(stack_fname{f,m}(end-7:end-4)) '...']);
                break;
            end
        end

        try
            for m=1 : size(stack_fname,2) 

                if ~isempty(stack_fname{f,m})
                    tic;
                    [refs{f,m}, num_frames(f,m)] = extract_candidate_reference_frames(stack_fname{f,m}, vertical_fringes_desinusoid_matrix', STRIP_SIZE, BAD_STRIP_THRESHOLD, MIN_NUM_FRAMES_PER_GROUP);
                    toc;
                end

            end
        catch ex
            warning on;
            warning(['Failed to find a reference frame in:' stack_fname{f,m}])
            warning(ex.message)
            warning(['From file: ' ex.stack(1).name ' Line: ' num2str(ex.stack(1).line)]);
        end

        %% Look for correspondence between all of the modalities.
        intersected = [];

        for m=1 : size(stack_fname,2)
            intersected = union(intersected, refs{f,m});
        end

        intersected(intersected==-1) = [];

        % First consider frames that appear in all 3 modalities.
        average_rank = nan(length(intersected),1);
        for r=1:length(intersected)
            of_interest = intersected(r);

            whichind = NaN(size(stack_fname,2) ,1);
            for m=1 : size(stack_fname,2)
                whichind(m) = NaN;% Negate a frame if it doesnt' appear in all modalities.
                rank = find( refs{f,m}==of_interest );
                if ~isempty(rank)
                    whichind(m) = rank*MODALITY_WEIGHTS(m);
                end
            end
            average_rank(r) = sum(whichind)./ sum(MODALITY_WEIGHTS);
        end

        [rankings, rankinds ] = sort(average_rank,1,'ascend');

        notallmodes_intersected = intersected(rankinds(isnan(rankings)));  
        intersected = intersected(rankinds(~isnan(rankings)));        
        
        % Then consider all the other frames.
        average_rank = nan(length(notallmodes_intersected),1);
        for r=1:length(notallmodes_intersected)
            of_interest = notallmodes_intersected(r);

            whichind = ones(size(stack_fname,2) ,1); 
            for m=1 : size(stack_fname,2)
                whichind(m) = length(refs{f,m})-length(intersected);% Weight heavily against a reference frame if it doesn't show in all modalities.
                rank = find( refs{f,m}==of_interest );
                if ~isempty(rank)
                    whichind(m) = rank*MODALITY_WEIGHTS(m);
                end
            end

            average_rank(r) = sum(whichind)./ sum(MODALITY_WEIGHTS);
        end

        [rankings, rankinds ] = sort(average_rank,1,'ascend');

        notallmodes_intersected = notallmodes_intersected(rankinds);
        
        intersected = [intersected; notallmodes_intersected];

        % Go through each intersected value and determine which group its in;
        % make separate rows in newrefs for disparate groups.
        grps = -ones(length(intersected),size(stack_fname,2));
        for i=1:length(intersected)
            for m=1 : size(stack_fname,2)
                [~, grp] = ind2sub( size(refs{f,m}), find(intersected(i)==refs{f,m}) );
                if ~isempty(grp)
                    grps(i,m) = grp;
                end
            end
        end

        newrefs = cell(1,100);    
        for m=1:size(grps,2)
            max_grp = max(grps(:,m));
            for g=1:max_grp

                ingrp = intersected( grps(:,m)==g );
                for n=1:length(newrefs)                
                    if isempty(newrefs{n})
                        newrefs{n} = ingrp;
                        break;
                    elseif ~isempty( intersect(newrefs{n}, ingrp ) )                    
                        newrefs{n} = [newrefs{n}; setdiff(ingrp,cell2mat(newrefs'))]; % If it already has been called into any other group, then don't include it in this one.
                        break;
                    end
                end

            end
        end

        newrefs = newrefs(~cellfun(@isempty,newrefs));

        for m=1:size(stack_fname,2)
            if ~isempty(stack_fname{f,m})
                vidnum = stack_fname{f,m}(end-7:end-4);
                break;
            end
        end

        %% Re-rank them based on their location in the intersected list.
        for g=1:length(newrefs)
            theserefs = newrefs{g};
            rankedrefs = -ones(size(theserefs));
            for i=1:length(theserefs)
                rankedrefs(i) = find(intersected==theserefs(i));
            end
            [rankings, rankinds ] = sort(rankedrefs);
            newrefs{g} = theserefs(rankinds);

            % Find out which reference frames fit best with which modalities
            if length(newrefs{g})>=NUM_REF_OUTPUT
                bestrefs = newrefs{g}(1:NUM_REF_OUTPUT)';
            else
                bestrefs = padarray(newrefs{g},[NUM_REF_OUTPUT-length(newrefs{g}) 0], NaN,'post')';
            end
            ref_best_modality = cell(size(bestrefs));
            ref_best_modality_inds = zeros(size(bestrefs));

            rel_ref = cell(size(refs(f,:)));% Compare each frame we've picked out relative to each other in their lists;
                                               % We only care if x is better than y, not if x is 3 indexes better than y
            for m=1:length(MODALITIES)
                if ~isempty(refs{f,m})
                    rel_ref_inds = find(sum(bestrefs==refs{f,m},2));
                    rel_ref{m} = refs{f,m}(rel_ref_inds);
                end
            end
            
            for r=1:length(bestrefs)
                thisrefrank = 100*ones(1,size(refs,2));
                for m=1:size(refs,2)
                    rank = find( rel_ref{m}==bestrefs(r) );
                    if ~isempty(rank)
                        thisrefrank(m) = rank*MODALITY_WEIGHTS(m);
                    end
                end
                [~, refrank_ind] = min(thisrefrank); % Whichever has the lowest index (best rank), record as the suggested modality.
                ref_best_modality{r} = MODALITIES{refrank_ind};
                ref_best_modality_inds(r) = refrank_ind;
            end

            % Write all of this to disk.
            fid= fopen(fullfile(mov_path{1},'Reference_Frames.csv'),'a');

            fprintf(fid,'"%s",',vidnum);

            for r=1:length(bestrefs)
                fprintf(fid,'"%s",%d,',ref_best_modality{r},bestrefs(r));
            end
            fprintf(fid,'\n');

            fclose(fid);

            for r=1:length(bestrefs)

                if ~isnan(bestrefs(r))
                    dmb_contents = default_dmb_contents;

                    dmb_contents.reference_frame = int32(bestrefs(r)-1);
                    for m=1:length(MODALITIES)
                        if m== ref_best_modality_inds(r)
                            [dmb_contents.image_sequence_absolute_path, dmb_contents.image_sequence_file_name]=getparent(stack_fname{f,m});
                            dmb_contents.n_frames = int32(num_frames(f,m));
                        elseif ~isempty(stack_fname{f,m})
                            [par, kid]=getparent(stack_fname{f,m});
                            dmb_contents.secondary_sequences_file_names = [dmb_contents.secondary_sequences_file_names; {kid}];
                            dmb_contents.secondary_sequences_absolute_paths = [dmb_contents.secondary_sequences_absolute_paths; {par}];                    
                        end
                    end
                    dmb_contents.user_defined_suffix = ['_ref_' num2str(bestrefs(r)) '_lps_' num2str(LPS) '_lbss_' num2str(LBSS) '_autogen' ];

                    if isempty(dmb_contents.secondary_sequences_file_names)
                       dmb_contents.secondary_sequences_file_names=''; 
                    end
                    if isempty(dmb_contents.secondary_sequences_absolute_paths)
                       dmb_contents.secondary_sequences_absolute_paths=''; 
                    end

    %                  save('test.mat','dmb_contents');
                    write_dmb_file(fullfile(mov_path{1}, [dmb_contents.image_sequence_file_name(1:end-4) dmb_contents.user_defined_suffix '.dmb']),dmb_contents);
                end
            end

        end


    end
end


       

close(h);




