function t3 = addterms(t1,t2,gradon)

    s1 = func2str(t1.fHandle);
    [vars,expr] = get_vars_str(s1);
    s2 = func2str(t2.fHandle);
    t3 = term('fHandle',eval([expr,strrep(s1,expr,''),'+',fcn_rep_vars(s2,vars)]),'gradon',0);

    t3.ftag = {t1.ftag,t2.ftag}; %%% need better way to combine these
    t3.nstates = t1.nstates;
    if gradon
        t3.gradterms = repmat(term('gradon',0),1,t3.nstates);
        for j=1:t3.nstates
            t3.gradterms(j) = addterms(t1.gradterms(j),t2.gradterms(j),0);
        end
    end

end