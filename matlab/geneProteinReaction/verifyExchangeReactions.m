function [exchangeRxns] = verifyExchangeReactions(model)
    exchangeRxns = model.rxns(findExcRxns(model));
end

