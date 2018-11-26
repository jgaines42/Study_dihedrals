%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_U_with_BSC
%
% Script takes unbound monomers and changes all side chains to the
% dihedral angles found in the bound structure
%
% Output:
%  _r_u_bSC.pdb : PDB file of unbound backbone with bound SC
%  _l_u_bSC.pdb
%
% Uses several functions from Single_rotation_code/ and PPI_analysis/
% Calls
% - fix_numbering_receptor_unbound
% - fix_numbering_ligand_unbound
% - swap_sidechains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/Users/jennifergaines/Documents/Code/Min_energy/Single_rotation_code/')
load('/Users/jennifergaines/Documents/Data/benchmark5/structures/PPI_names_new.mat')
addpath('/Users/jennifergaines/Documents/Code/PPI_analysis/')

% Remove AB structures since don't have 2 unique structures
ind0 = ismember(names(:,2),{'AB'});
to_use = ~ind0;
names(~ind0,8) = num2cell(1);


dih = [];
%tot_names = {};
%mutations = [];
%start = i;
%missing = {};
%cant_run = {};
%1:84,86:154,156, 158:160,162:177,
for i =[1:size(names,1)] % Skipped 85, 155, 157, 161, 168,  178, 218
    i
    if to_use(i) == 1
        
        %Load combined bound structure
        load(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/',names{i,1}, '_combined.mat'))
        orig = tempModel2;
        
        % Figure out numbering for complex and interface data relative to raw data
        [~,x] = unique(tempModel2(:,5));
        chains = tempModel2(sort(x),5); %Unique chains, in order found in PDB
        res_ids = [];
        all_ids = cell2mat(tempModel2(:,6));
        chain_protein = cell(size(chains,1),1);
        ids = cell(size(chains,1),1);
        
        %Loop over posible chains
        for j = 1:size(chains,1)
            ind0 = ismember(tempModel2(:,5), chains{j});
            ind1 = ismember(tempModel2(:,2), {'CA'});
            res_1 = cell2mat(tempModel2(ind0&ind1,6));
            res_ids = [res_ids; res_1,repmat(j,size(res_1,1),1)];
            all_ids(ind0,2) = j;
            chain_protein(j,1) = {tempModel2(ind0,:)};
            ids(j) = {cell2mat(tempModel2(ind0,6))};
            
        end
        % ids: cell array, one cell for each chain, containing all column 6 (ids) of tempModel2 for each chain
        % chain_protein: one cell for each chain, containing all of tempModel2
        % all_ids: column 1: residue ids, column 2: index of chain (1, 2, etc)
        % res_ids: one row for each residues, column 1: res id, column 2: index of chain
        
        complex = tempModel2;
        
        %% Run for unbound monomers
        if exist(strcat('/Users/jennifergaines/Documents/Code/Alter_sidechains/', names{i,1}, '_r_u_fixed_H.pdb'), 'file')
            pdbstruct = pdbread(strcat('/Users/jennifergaines/Documents/Code/Alter_sidechains/', names{i,1}, '_r_u_fixed_H.pdb'));
            x = pdbstruct.Model.Atom;
            tempModel1=struct2cell(x);
            tempModel2=reshape(tempModel1,size(tempModel1,1),size(tempModel1,3))';
            save(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/monomers/', names{i,1}, '_r_H_fixed.mat'), 'tempModel2')
            load(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/monomers/', names{i,1}, '_r_H_fixed.mat'))  
        else
            load(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/monomers/', names{i,1}, '_r_u_H.mat'))
        end
       %{
       if ismember(i, [6,9, 20,26, 36,47, 52, 54,58,61, 62,63,77, 79, 82, 85,91, 100, 108, 109, 112, 116 122 123 129 132 133])
        end
       %}
        
        orig_receptor = tempModel2;
         [~,x] = unique(tempModel2(:,5));
       orig_chain = tempModel2(sort(x),5);
        %Fix numbering if needed to be the same as bound complex
        tempModel2 = fix_numbering_receptor_unbound(i, tempModel2, orig);
        number_mapping = [cell2mat(orig_receptor(:,6)), cell2mat(tempModel2(:,6))];
        %Get chains
        [~,x] = unique(tempModel2(:,5));
        c1 = tempModel2(sort(x),5);
        res_ids_r = [];
        all_ids_r = cell2mat(tempModel2(:,6));
        ids = cell(size(c1,1),1);
        new_chain = cell(size(tempModel2,1),1);
        for j = 1:size(c1,1)
            ind0 = ismember(tempModel2(:,5), c1{j});
            res_1 = unique(cell2mat(tempModel2(ind0,6)));
            res_ids_r = [res_ids_r; res_1,repmat(j,size(res_1,1),1)];
            new_chain(ind0) = chains(j);
            all_ids_r(ind0,2) = j;
            x = tempModel2(ind0,:);
            x(:,5) = chains(j);
            chain_protein(j,2) = {x};
            ids(j,1) = {cell2mat(tempModel2(ind0,6))};
            ids(j,2) = {number_mapping(ind0,:)};
        end
        tempModel2(:,5) = new_chain;
        receptor = tempModel2;
        c1 = unique(receptor(:,5));
        
        start_chain = j+1;
        
        
        %% Preprocess ligand
        
        if exist(strcat('/Users/jennifergaines/Documents/Code/Alter_sidechains/', names{i,1}, '_l_u_fixed_H.pdb'), 'file')
            pdbstruct = pdbread(strcat('/Users/jennifergaines/Documents/Code/Alter_sidechains/', names{i,1}, '_l_u_fixed_H.pdb'));
            x = pdbstruct.Model.Atom;
            tempModel1=struct2cell(x);
            tempModel2=reshape(tempModel1,size(tempModel1,1),size(tempModel1,3))';
            save(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/monomers/', names{i,1}, '_l_H_fixed.mat'), 'tempModel2')
            load(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/monomers/', names{i,1}, '_l_H_fixed.mat'))  
        else  
                    load(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/monomers/', names{i,1}, '_l_u_H.mat'))
        end

        orig_ligand = tempModel2;
          [~,x] = unique(tempModel2(:,5));
       orig_chain_l = tempModel2(sort(x),5);
        %Fix numbering if needed to be the same as bound complex
        [tempModel2] = fix_numbering_ligand_unbound(i, tempModel2);
        if ismember(names{i,1}, '1JMO') %1JMO gets chains swapped
            ind0 = ismember(orig_ligand(:,5),'L');
            orig_ligand = [orig_ligand(ind0,:);orig_ligand(~ind0,:)];
        end
        number_mapping = [cell2mat(orig_ligand(:,6)), cell2mat(tempModel2(:,6))];

        if ismember(names{i,1}, '1FC2')
            ind0 = ismember(tempModel2(:,5), 'A');
            tempModel2 = tempModel2(ind0,:);
        elseif ismember(names{i,1}, '2VIS')
            ind0 = ismember(tempModel2(:,5), 'A');
            tempModel2 = tempModel2(ind0,:);
        end
        
        %Get chain data
        [~,x] = unique(tempModel2(:,5));
        c2 = tempModel2(sort(x),5);
        res_ids_l = [];
        all_ids_l = cell2mat(tempModel2(:,6));
        new_chain = cell(size(tempModel2,1),1);
        for j = 1:size(c2,1)
            ind0 = ismember(tempModel2(:,5), c2{j});
            res_1 = unique(cell2mat(tempModel2(ind0,6)));
            res_ids_l = [res_ids_l; res_1,repmat(start_chain,size(res_1,1),1)];
            new_chain(ind0,1) = chains(start_chain);
            all_ids_l(ind0,2) = start_chain;
            x = tempModel2(ind0,:);
            x(:,5) = chains(start_chain);
            chain_protein(start_chain,2) = {x};
            ids(start_chain,1) = {cell2mat(tempModel2(ind0,6))};
            ids(start_chain,2) = {number_mapping(ind0,:)};

            start_chain = start_chain+1;

        end
        tempModel2(:,5) = new_chain;
        ligand = tempModel2;
        
        %Combine receptor and ligand
        tempModel2 = [receptor;ligand];
        % Create unique index in first column
        
        tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
        %Separate out to just have one row per AA
        ind0 = ismember(tempModel2(:,2) , 'CA');
        sub_protein = tempModel2(ind0,:);
        sub_ids = cell2mat(sub_protein(:,6));
        sub_chains = sub_protein(:,5);
        
        new_protein = {};
        orig_protein = {};
        
        % Loop over all chains
        for j = 1:size(chains,1)
            unbound_chain = chain_protein{j,2};
            unbound_chain(:,1) = num2cell([1:size(unbound_chain,1)]);
            bound_chain = chain_protein{j,1};
            bound_chain(:,1) = num2cell([1:size(bound_chain,1)]);
            this_mapping = ids{j,2};
            this_id_u = cell2mat(unbound_chain(:,6));
            this_id_b = cell2mat(bound_chain(:,6));
            res_ids = unique(this_id_u);
            
            % Loop over all residue in unbound structure
            for residues = 1:size(res_ids,1)
                %Set up the residue
                ind0 = find(this_id_u ==res_ids(residues));
                resiName = unbound_chain{ind0(1),4}; %Get name of amino acid
                resiName(2:3) = lower(resiName(2:3));
                
                % If this residue exists in the bound structure
                if sum(ismember(this_id_b, res_ids(residues)))>=1
                    
                    %Isolate dipeptide of each structure
                    [allDipeptide_u, next_pro, numAtom, DOF] = isolate_dipeptide(unbound_chain, this_id_u, res_ids(residues), resiName);
                    [allDipeptide_b, next_pro1, numAtom1, DOF1] = isolate_dipeptide(bound_chain, this_id_b, res_ids(residues), resiName);
                    

                    % If they are both are *not* terminal residues
                    if sum(ismember(allDipeptide_u(:,2), {'OXT'})) == 0 &  sum(ismember(allDipeptide_b(:,2), {'OXT'}))==0
                       
                        % If there is no N-1 residue in the bound structure, pad with blank cells
                        if cell2mat(allDipeptide_b(1,6)) == res_ids(residues)
                            allDipeptide_b = [cell(3,16);allDipeptide_b];
                            allDipeptide_b(1:3,14) = {'O'};
                            allDipeptide_b(1:3,8:10) =num2cell(0);
                            allDipeptide_b(1:3,6) = {res_ids(residues)-1};
                        end
                        
                        % If there is no N-1 residue in the unbound structure, pad with blank cells
                        if cell2mat(allDipeptide_u(1,6)) == res_ids(residues)
                            allDipeptide_u = [cell(3,16);allDipeptide_u];
                            allDipeptide_u(1:3,14) = {'O'};
                            allDipeptide_u(1:3,8:10) =num2cell(0);
                            allDipeptide_u(1:3,6) = {res_ids(residues)-1};
                        end
                        
                        % If there is no N+1 residue in the unbound structure, pad with blank cells 
                        if cell2mat(allDipeptide_u(size(allDipeptide_u,1),6)) == res_ids(residues)
                            s = size(allDipeptide_u,1);
                            allDipeptide_u = [allDipeptide_u;cell(3,16)];
                            allDipeptide_u(s+1:s+2,14) = {'O'};
                            allDipeptide_u(s+3,14) = {'H'};
                            allDipeptide_u(s+1:s+3,8:10) =num2cell(0);
                            allDipeptide_u(s+1:s+3,6) = {res_ids(residues)+1};
                            next_pro = 0;
                        end
                        
                        % If there is no N+1 residue in the bound structure, pad with blank cells 
                        if cell2mat(allDipeptide_b(size(allDipeptide_b,1),6)) == res_ids(residues)
                            s = size(allDipeptide_b,1);
                            allDipeptide_b = [allDipeptide_b;cell(3,16)];
                            allDipeptide_b(s+1:s+2,14) = {'O'};
                            allDipeptide_b(s+3,14) = {'H'};
                            allDipeptide_b(s+1:s+3,8:10) =num2cell(0);
                            allDipeptide_b(s+1:s+3,6) = {res_ids(residues)+1};
                            next_pro1 = 0;
                        end
                        
                        % Check dipeptide order of unbound structure
                        if next_pro == 0
                            [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide_u(4:size(allDipeptide_u,1)-3,:), resiName);
                            allDipeptide_u(4:size(allDipeptide_u,1)-3,:)= new_Dipeptide;
                        else
                            [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide_u(4:size(allDipeptide_u,1)-2,:), resiName);
                            allDipeptide_u(4:size(allDipeptide_u,1)-2,:)= new_Dipeptide;
                        end
                        
                        % Check dipeptide order of bound structure
                        if next_pro1 == 0
                            [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide_b(4:size(allDipeptide_b,1)-3,:), resiName);
                            allDipeptide_b(4:size(allDipeptide_b,1)-3,:)= new_Dipeptide;
                        else
                            [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide_b(4:size(allDipeptide_b,1)-2,:), resiName);
                            allDipeptide_b(4:size(allDipeptide_b,1)-2,:)= new_Dipeptide;
                        end
                        
                        % remove H
                        ind0 = ismember(allDipeptide_b(:,14), 'H');
                        allDipeptide_b = allDipeptide_b(~ind0,:);
                        ind0 = ismember(allDipeptide_u(:,14), 'H');
                        allDipeptide_u = allDipeptide_u(~ind0,:);
                        
                        % If dipeptides are now the same size, with the same residue name, swap the side chains
                        if size(allDipeptide_b,1) == size(allDipeptide_u,1) && strcmp(allDipeptide_b(4,4),allDipeptide_u(4,4))
                            [new_dipeptide] = swap_sidechain(allDipeptide_u, allDipeptide_b,resiName, next_pro);
                            allDipeptide_u = new_dipeptide;
                            ind0 = ismember(cell2mat(allDipeptide_u(:,6)), res_ids(residues));
                            new_protein = [new_protein;allDipeptide_u(ind0,:)];
                            ind1 = ismember(cell2mat(allDipeptide_b(:,6)), res_ids(residues));
                            orig_protein = [orig_protein;allDipeptide_b(ind1,:)];
                        else 
                            % See if there is a mutation
                            if ~strcmp(allDipeptide_b(4,4),allDipeptide_u(4,4))
                                mutations = [mutations;i,res_ids(residues),j];
                            % If no mutation, must have missing atoms
                            % Save to 'missing' to alter using pymol
                            elseif size(allDipeptide_u,1) <= size(allDipeptide_b,1)
                               ind0 = ismember(this_mapping(:,2), res_ids(residues));
                               map_1 = this_mapping(ind0,1);
                               map_2 = unique(map_1);
                               if size(map_2,1) > 1
                                   here = 1;
                               end
                                if j <= size(c1,1)
                                    missing = [missing; names(i,1), num2cell(map_2),orig_chain(j),'r' , allDipeptide_b(5,4)];
                                else
                                    missing = [missing; names(i,1), num2cell(map_2),orig_chain_l(j-size(c1,1)),'l' , allDipeptide_b(5,4)];
                                end
                            end
                            % Because the two dipeptides aren't the same,
                            % just keep original side chain
                            ind0 = ismember(cell2mat(allDipeptide_u(:,6)), res_ids(residues));
                            new_protein = [new_protein;allDipeptide_u(ind0,:)];
                            ind1 = ismember(cell2mat(allDipeptide_b(:,6)), res_ids(residues));
                            orig_protein = [orig_protein;allDipeptide_b(ind1,:)];
                        end
                    else % If this was the end of the protein chain, keep the same side chain
                        here = 1;
                        ind0 = ismember(cell2mat(allDipeptide_u(:,6)), res_ids(residues));
                        new_protein = [new_protein;allDipeptide_u(ind0,:)];
                        ind1 = ismember(cell2mat(allDipeptide_b(:,6)), res_ids(residues));
                        orig_protein = [orig_protein;allDipeptide_b(ind1,:)];
                        
                    end
                end
            end
        end
        
        ind0 = ismember(new_protein(:,5), c1);
        new_receptor = new_protein(ind0,:);
        new_ligand = new_protein(~ind0,:);
        
        if size(new_receptor,1) == 0 || size(new_ligand,1) == 0
           display('problem with %s', names{i,1})
        end
        
        % Save to files
        f = fopen(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/modified_SC/',names{i,1}, '_r_u_bSC.pdb'), 'w');
        tempModel2 = new_receptor;
        for atoms = 1:size(tempModel2,1)
            if size(tempModel2{atoms,2},2)<4
                fprintf(f,( '%-6s%5d  %-3s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f         \n'), 'ATOM', atoms, tempModel2{atoms,2}, tempModel2{atoms,4}, tempModel2{atoms,5}, tempModel2{atoms,6}, '', tempModel2{atoms,8}, tempModel2{atoms,9}, tempModel2{atoms,10}, tempModel2{atoms,11}, tempModel2{atoms,11});
            else
                fprintf(f, ('%-6s%5d %-4s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f          \n'), 'ATOM', atoms, tempModel2{atoms,2}, tempModel2{atoms,4},tempModel2{atoms,5}, tempModel2{atoms,6}, '', tempModel2{atoms,8}, tempModel2{atoms,9}, tempModel2{atoms,10}, tempModel2{atoms,11}, tempModel2{atoms,11});
            end
        end
        fclose(f);
        f = fopen(strcat('/Users/jennifergaines/Documents/Data/benchmark5/structures/modified_SC/',names{i,1}, '_l_u_bSC.pdb'), 'w');
        tempModel2 = new_ligand;
        for atoms = 1:size(tempModel2,1)
            if size(tempModel2{atoms,2},2)<4
                fprintf(f,( '%-6s%5d  %-3s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f         \n'), 'ATOM', atoms, tempModel2{atoms,2}, tempModel2{atoms,4}, tempModel2{atoms,5}, tempModel2{atoms,6}, '', tempModel2{atoms,8}, tempModel2{atoms,9}, tempModel2{atoms,10}, tempModel2{atoms,11}, tempModel2{atoms,11});
            else
                fprintf(f, ('%-6s%5d %-4s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f          \n'), 'ATOM', atoms, tempModel2{atoms,2}, tempModel2{atoms,4},tempModel2{atoms,5}, tempModel2{atoms,6}, '', tempModel2{atoms,8}, tempModel2{atoms,9}, tempModel2{atoms,10}, tempModel2{atoms,11}, tempModel2{atoms,11});
            end
        end
        fclose(f);
        
        
    end
end
