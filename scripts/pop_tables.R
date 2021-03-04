#### SETTINGS ####
####################################################################################################!

#### MAIN ####
####################################################################################################!

## -- Dataset S1 -- ##

ds1_pop_df <- readxl::read_xls("dataset_S1/data/populations/populations.xls") %>% as.data.table()

ds1_natives <- ds1_pop_df[Selection_Group %in% c("Amazonia", "Mesoamerica"), .(Population, Macro_language, Region, Latitude, Longitude, Source, N)] %>%
  rename(MacroLanguage = "Macro_language") %>%
  mutate_if(is.numeric, round, 2) %>%
  arrange(Region, Population, MacroLanguage, Source) %>%
  # Translate to Portuguese
  set_names(c("População", "Grupo Linguístico", "Região",
              "Latitude", "Longitude", "Referência", "N")) %>%
  mutate(Referência=gsub("This", "Este", Referência),
         Referência=gsub("Castro_Silva", "Castro e Silva", Referência),
         Região=case_when(Região=="Brazil"~"Brasil", Região=="Mexico"~"México", T~"NA"))

ds1_selgroups <- ds1_pop_df[!is.na(Selection_Group), .(N = sum(N), Macroregion = unique(Macro_region)), by = Selection_Group] %>%
  rename(Grupo=Selection_Group, Macrorregião=Macroregion) %>%
  select(Grupo, Macrorregião, N) %>%
  mutate (
    Grupo=case_when(Grupo=="Amazonia"~"Amazônia", Grupo=="EastAsia"~"Leste da Ásia", Grupo=="Mesoamerica"~"Mesoamérica", T~"NA"),
    Macrorregião=case_when(Macrorregião=="South_America"~"América do Sul", Macrorregião=="East_Asia"~"Ásia", Macrorregião=="Central_America"~"América Central", T~"NA")
  )

ds1_natives_kbl <- kbl(ds1_natives, booktabs = T, format = "latex", linesep = "", escape = T) %>%
  kable_styling(latex_options = c("striped", "scale_down"))

ds1_selgroups_kbl <- kbl(ds1_selgroups, booktabs = T, format = "latex", linesep = "", escape = F) %>%
  kable_classic()

save_kable(ds1_natives_kbl, "tables/tex/ds1_natives.tex")
save_kable(ds1_selgroups_kbl, "tables/tex/ds1_selgroups.tex")

## -- Dataset S2 -- ##

ds2_pop_df <- fread("dataset_S2/data/populations/populations.csv")
ds2_pop_df[grepl("Yanesha", FID), FID := "Yanesha"]

ds2_pop_df[sel_group=='Amazonia', REGION := "South_America"]
ds2_pop_df[sel_group=='Andes', REGION := "South_America"]
ds2_pop_df[sel_group=='Mesoamerica', REGION := "Mesoamerica"]
ds2_pop_df[sel_group=='EastAsia', REGION := "East_Asia"]

ds2_natives <- ds2_pop_df %>%
  select(FID, COUNTRY, sel_group, REGION, SOURCE) %>%
  filter(sel_group %in% c("Amazonia", "Mesoamerica")) %>%
  add_count(FID, name = "N") %>%
  distinct(FID, .keep_all = TRUE) %>%
  arrange(sel_group, FID, COUNTRY, -N) %>%
  rename(População=FID, País=COUNTRY, Região=REGION, Referência=SOURCE, Ecorregião=sel_group) %>%
  mutate (
    País=case_when(País=="Brazil"~"Brasil", País=="Mexico"~"México", T~"Peru"),
    Ecorregião=case_when(Ecorregião=="Amazonia"~"Amazônia", Ecorregião=="Mesoamerica"~"Mesoamérica", T~"NA"),
    Região=case_when(Região=="South_America"~"América do Sul", Região=="Mesoamerica"~"América Central", T~"NA")
  )

ds2_selgroups <- ds2_pop_df[sel_group != "", .(N = .N, Macroregion = unique(REGION)), by = sel_group] %>%
  filter(sel_group!="Andes") %>%
  rename(Grupo=sel_group, Macrorregião=Macroregion) %>%
  select(Grupo, Macrorregião, N) %>%
  mutate (
    Grupo=case_when(Grupo=="Amazonia"~"Amazônia", Grupo=="EastAsia"~"Leste da Ásia", Grupo=="Mesoamerica"~"Mesoamérica", T~"NA"),
    Macrorregião=case_when(Macrorregião=="South_America"~"América do Sul", Macrorregião=="East_Asia"~"Ásia", Macrorregião=="Mesoamerica"~"América Central", T~"NA")
  )

ds2_natives_kbl <- kbl(select(ds2_natives, -Região), booktabs = T, format = "latex", linesep = "", escape = F) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  pack_rows("América do Sul", 1, 7) %>%
  pack_rows("Mesoamérica", 8, nrow(ds2_natives))

ds2_selgroups_kbl <- kbl(ds2_selgroups, booktabs = T, format = "latex", linesep = "", escape = F) %>%
  kable_classic()

save_kable(ds2_natives_kbl, "tables/tex/ds2_natives.tex")
save_kable(ds2_selgroups_kbl, "tables/tex/ds2_selgroups.tex")

# XLS

xls <- list (
  ds1_natives=ds1_natives, ds1_selgroups=ds1_selgroups,
  ds2_natives=ds2_natives, ds2_selgroups=ds2_selgroups
  )

WriteXLS::WriteXLS(xls, "tables/xlsx/populations.xls", AdjWidth = T, BoldHeaderRow = T)
