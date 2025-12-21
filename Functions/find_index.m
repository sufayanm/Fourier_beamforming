% This function finds the index of closest number to 'dist' in 'in_axis'

function out_index = find_index(dist, in_axis, out_range)

if nargin<3
    out_range = 0;
end

if isscalar(dist)

    if dist>max(in_axis)||dist<min(in_axis)
        error("Querry value is not in the range of sample points")
    end
    [~, out_index] = min(abs(dist - in_axis(:))) ;
else
    size_dist = size(dist) ;
    dist = reshape(dist, 1, []) ;
    in_axis = in_axis(:) ;

    if max(dist)>max(in_axis)||min(dist)<min(in_axis)
        if out_range==0
            error("Solution :- Some of the querry values are not in the range of sample points. " + ...
                "If you want to proceed anyway use 1 as a third argument while calling the function")
        end
    end
    [~, out_index] = min(abs(dist - in_axis)) ;
    out_index =  out_index(:) ;

    out_index = reshape(out_index, size_dist) ;
end



