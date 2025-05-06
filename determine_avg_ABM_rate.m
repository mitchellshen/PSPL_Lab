function [avg_IN,avg_OUT,avg_COIN,avg_ADE,avg_IN_err,avg_OUT_err,avg_COIN_err,avg_ADE_err,t_desired,t_elapsed] = determine_avg_ABM_rate( ...
    tdesired1,tdesired2,t,IN_CEM_rate,OUT_CEM_rate,COIN_rate)
    indt = t>=tdesired1 & t <=tdesired2;
    t_desired = t(indt); 
    t_elapsed = seconds(t_desired(end) - t_desired(1));
    avg_IN = sum(IN_CEM_rate(indt))/t_elapsed;
    avg_IN_err = sqrt(avg_IN/t_elapsed);
    avg_OUT = sum(OUT_CEM_rate(indt))/t_elapsed;
    avg_OUT_err = sqrt(avg_OUT/t_elapsed);
    avg_COIN = sum(COIN_rate(indt))/t_elapsed;
    avg_COIN_err = sqrt(avg_COIN/t_elapsed);
    avg_ADE = (avg_COIN.^2)./(avg_IN.*avg_OUT);
    avg_ADE_err = sqrt((2*avg_COIN*avg_COIN_err/(avg_IN.*avg_OUT)).^2 + (avg_COIN.^2/(avg_IN.^2*avg_OUT)*avg_IN_err).^2 + ...
    + (avg_COIN.^2/(avg_OUT.^2*avg_IN)*avg_OUT_err).^2);
end