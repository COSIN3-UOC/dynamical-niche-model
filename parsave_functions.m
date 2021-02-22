classdef parsave_functions
%{
To be able to save files in in the main file when running realization loops
with parfor
%}  
    methods(Static)
          function parsave_matrices( fname, out_matrix)
            save( fname, 'out_matrix');
          end
          
          function parsave_niche_overlaps( fname, niche_matrix)
            save( fname, 'niche_matrix');
          end
          
          function parsave_abundances( fname, tTot, yTot)
            save( fname, 'tTot', 'yTot');
          end
          
          function parsave_niche_centers( fname, Ha, Hp)
            save( fname, 'Ha', 'Hp');
          end

    end
end