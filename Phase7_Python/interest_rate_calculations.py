Initial_principal =100000
interest_rate = 0.055
payment_rate = 0.1

years_of_investment = 7
principal = Initial_principal
cash_withdrawal = 0
for i in range(years_of_investment-1):
    difference = round(principal * (1 + interest_rate) - principal,2)
    principal = round(principal * (1 + interest_rate),2) 
    cash_withdrawal_this_round = round(principal * payment_rate,2)
    cash_withdrawal += cash_withdrawal_this_round
    print("End of Year ",i+1,", total: ", principal, ", taxable_interest", difference,", Cash withdrawal: " ,cash_withdrawal_this_round,", Principal for next year: ", round(principal - principal * payment_rate,2))
    principal = round(principal - principal * payment_rate,2)

difference = round(principal * (1 + interest_rate) - principal,2)
principal = round(principal * (1 + interest_rate),2)
print("End of Year: ", 7, "Total: ", principal, ", taxable_interest", difference)
print("Cash withdrawal: ", round(cash_withdrawal,2))
print("Final principal: ", principal)
print("Total: ", cash_withdrawal + principal)

print("Total if invested in ",years_of_investment," years and only take out interest: ", round(Initial_principal * (1 + interest_rate)**years_of_investment,2))


print("Total if invested in ",years_of_investment," years without taking any money out: ", round(Initial_principal * (1 + interest_rate)**years_of_investment,2))