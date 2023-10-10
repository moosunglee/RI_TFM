function Ref = Pad_tomograms(Ref,dims,RI_bg)

        if ~isempty(Ref.I_WF)
            for j1 =1:3
                if ~isempty(Ref.I_SIM{j1})
                    Ref.I_SIM{j1} = My_paddzero_MS(Ref.I_SIM{j1},dims);
                end
                if ~isempty(Ref.I_WF{j1})
                    Ref.I_WF{j1} = My_paddzero_MS(Ref.I_WF{j1},dims);
                end
            end
        end
            if isfield(Ref, 'RI_multi')
                Ref.RI_multi = My_paddzero_MS(Ref.RI_multi,dims, RI_bg);
            end
            if isfield(Ref, 'RI_stitched')
                Ref.RI_stitched = My_paddzero_MS(Ref.RI_stitched,dims, RI_bg);
            end

end