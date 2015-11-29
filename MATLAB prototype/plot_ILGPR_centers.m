function plot_ILGPR_centers(ilgpr,Interpolant,figureHandle)

    figure(figureHandle);

    for j = 1:ilgpr.nLGPs    
        if ilgpr.LGPs{j}.started == 1
            markerString = 'ms';
        else
            markerString = 'mx';
        end
        plot3(ilgpr.LGPs{j}.u(1),ilgpr.LGPs{j}.u(2),Interpolant(ilgpr.LGPs{j}.u(1),ilgpr.LGPs{j}.u(2)),markerString,'MarkerSize',10,'LineWidth',5);
    end

end