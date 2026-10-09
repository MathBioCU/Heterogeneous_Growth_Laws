



Dataset = {"ChickenData.csv","LoblollyData.csv","OrangeTreeData.csv"};
SaveName = {"ChickenData","LoblollyData","OrangeTreeData"};
Choice = 1;





Data = convert_nlme_csv_to_wsindy(Dataset{Choice});

save(SaveName{Choice},"Data")