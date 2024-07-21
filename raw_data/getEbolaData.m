function [initial_date, incidence_data] = getEbolaData()

incidence_data = zeros(100,1);

initial_date = '2-apr-2018'; % This date was chosen so that, when the data are aggregated into weeks, we get the (publically available) WHO weekly incidence curve

first_case = '5-apr-2018';

incidence_data(daysact(initial_date,  '5-apr-2018') + 1) = incidence_data(daysact(initial_date,  '5-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '8-apr-2018') + 1) = incidence_data(daysact(initial_date,  '8-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '12-apr-2018') + 1) = incidence_data(daysact(initial_date,  '12-apr-2018') + 1) + 3;
incidence_data(daysact(initial_date,  '13-apr-2018') + 1) = incidence_data(daysact(initial_date,  '13-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '18-apr-2018') + 1) = incidence_data(daysact(initial_date,  '18-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '19-apr-2018') + 1) = incidence_data(daysact(initial_date,  '19-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '20-apr-2018') + 1) = incidence_data(daysact(initial_date,  '20-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '21-apr-2018') + 1) = incidence_data(daysact(initial_date,  '21-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '23-apr-2018') + 1) = incidence_data(daysact(initial_date,  '23-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '24-apr-2018') + 1) = incidence_data(daysact(initial_date,  '24-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '25-apr-2018') + 1) = incidence_data(daysact(initial_date,  '25-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '27-apr-2018') + 1) = incidence_data(daysact(initial_date,  '27-apr-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '1-may-2018') + 1) = incidence_data(daysact(initial_date,  '1-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '3-may-2018') + 1) = incidence_data(daysact(initial_date,  '3-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '4-may-2018') + 1) = incidence_data(daysact(initial_date,  '4-may-2018') + 1) + 6;
incidence_data(daysact(initial_date,  '5-may-2018') + 1) = incidence_data(daysact(initial_date,  '5-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '6-may-2018') + 1) = incidence_data(daysact(initial_date,  '6-may-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '7-may-2018') + 1) = incidence_data(daysact(initial_date,  '7-may-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '8-may-2018') + 1) = incidence_data(daysact(initial_date,  '8-may-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '10-may-2018') + 1) = incidence_data(daysact(initial_date,  '10-may-2018') + 1) + 3;
incidence_data(daysact(initial_date,  '12-may-2018') + 1) = incidence_data(daysact(initial_date,  '12-may-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '13-may-2018') + 1) = incidence_data(daysact(initial_date,  '13-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '14-may-2018') + 1) = incidence_data(daysact(initial_date,  '14-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '15-may-2018') + 1) = incidence_data(daysact(initial_date,  '15-may-2018') + 1) + 3;
incidence_data(daysact(initial_date,  '16-may-2018') + 1) = incidence_data(daysact(initial_date,  '16-may-2018') + 1) + 3;
incidence_data(daysact(initial_date,  '18-may-2018') + 1) = incidence_data(daysact(initial_date,  '18-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '19-may-2018') + 1) = incidence_data(daysact(initial_date,  '19-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '20-may-2018') + 1) = incidence_data(daysact(initial_date,  '20-may-2018') + 1) + 3;
incidence_data(daysact(initial_date,  '21-may-2018') + 1) = incidence_data(daysact(initial_date,  '21-may-2018') + 1) + 2;
incidence_data(daysact(initial_date,  '28-may-2018') + 1) = incidence_data(daysact(initial_date,  '28-may-2018') + 1) + 1;
incidence_data(daysact(initial_date,  '2-jun-2018') + 1) = incidence_data(daysact(initial_date,  '2-jun-2018') + 1) + 1;