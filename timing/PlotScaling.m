function tf = PlotScaling(err_v_time,err_prob)

% function PlotScaling(err_v_time, err_prob) determines the minimum  
% time required to make the error < err_prob and plots tf vs err_prob 
% in a loglog plot
%
% input: if err_vs_time must be a structure with at least two fields
%        corresponding to the time (t) and error probability.
%        one field must be t, the name of the other fields are arbitrary.
% 
% output: tf (vector) or structure with same fields as err_v_time
%         containing the minimum times

if ~isfield(err_v_time,'t')
    error('time vector must be specified')
end
% get all fields ~= t
FIELDS = fields(err_v_time);
FIELDS = {FIELDS{cellfun(@(x)~strcmp(x,'t'),FIELDS)}};

for k=1:length(FIELDS)
    tmp = arrayfun(@(Err)err_v_time.t(find(err_v_time.(FIELDS{k})<Err,1)),err_prob,'UniformOutput',false);
    tf.(FIELDS{k}) = cell2mat(cellfun(@(x) x(~isempty(x)), tmp,'UniformOutput',false));
    n = length(tf.(FIELDS{k}));

    Fit.(FIELDS{k}) = polyfit(log10(err_prob(1:n)'),log10(tf.(FIELDS{k})'),1)
    
    figure(2)
    loglog(err_prob(1:n),tf.(FIELDS{k}))
    grid on, hold on
    xlabel('error level \epsilon_{prob} (%)')
    ylabel('minimum time required (1/J)')
end
% optional plot fit
for k=1:length(FIELDS)
    semilogy(err_prob,10.^polyval(Fit.(FIELDS{k}),log10(err_prob)))
    text(1e-4,1,sprintf('t_f = %0.2f \\epsilon_{prob}^{%0.2f}', 10^Fit.(FIELDS{k})(2),Fit.(FIELDS{k})(1)),'FontSize',14)
end

legend(FIELDS)
