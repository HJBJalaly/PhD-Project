function SavePdfFast(Name)

% saveas(gcf,[Name,'.fig'])
set(gcf,'Units','centimeters','PaperUnits','centimeters');
Pos=get(gcf,'Position');
set(gcf,'PaperPosition', [0 0 Pos(3) Pos(4)],'PaperSize', [Pos(3) Pos(4)])
% set(gcf,'PaperUnits',Units);
saveas(gcf,[Name,'.pdf'])
% system(['/usr/local/texlive/2014/bin/x86_64-linux/pdfcrop ',Name,'.pdf ',Name,'.pdf']);