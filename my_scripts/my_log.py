#! /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_26/external/slc7_amd64_gcc820/bin/python3

import pandas as pd

print()
df = pd.read_csv("~/CMSSW_10_6_26/src/tttt/test/output/all_data/2018/all_data/first_run/output_csv.csv")

df_os, df_ss = df[df["has_OS"] == 1], df[df["has_OS"] == 0]
os_events, ss_events = df_os.shape[0], df_ss.shape[0]
print(f"The total number of OS events: {os_events}")
print(f"The total number of SS events: {ss_events}")

assert os_events + ss_events == df.shape[0]

total_ss_gen = df["has_SS_gen"].value_counts().loc[1]
print(f"Total number of SS gen events: {total_ss_gen}")

os_ss_gen,ss_ss_gen = df_os["has_SS_gen"].value_counts().loc[1], df_ss["has_SS_gen"].value_counts().loc[1]
print(f"OS reco events with SS gen: {os_ss_gen}")
print(f"SS reco events with SS gen: {ss_ss_gen}")

assert os_ss_gen + ss_ss_gen == df["has_SS_gen"].value_counts().loc[1]

total_fakes = df["has_fake"].value_counts().loc[1]
print(f"The total number events with fakes: {total_fakes}")

os_fakes, ss_fakes = df_os["has_fake"].value_counts().loc[1], df_ss["has_fake"].value_counts().loc[1]
print(f"OS reco events with fakes: {os_fakes}") 
print(f"SS reco events with fakes: {ss_fakes}")

assert os_fakes + ss_fakes == total_fakes

df_two_genmatched = df_os[df_os["nGenMatched_leptons"] == 2]
print("For OS reco events")
for i in sorted(df_two_genmatched["nFlips"].value_counts().index):
	print(f"\tNumber of events with {i} flips: {df_two_genmatched['nFlips'].value_counts().loc[i]}")

df_two_genmatched = df_ss[df_ss["nGenMatched_leptons"] == 2]
print("For SS reco events")
for i in sorted(df_two_genmatched["nFlips"].value_counts().index):
	print(f"\tNumber of events with {i} flips: {df_two_genmatched['nFlips'].value_counts().loc[i]}")

df_eta = df[(df["leading_eta"] > 2.5) | (df["trailing_eta"] > 2.5)]
print(f"Events with eta greater than 2.5: {df_eta.shape[0]}")

max_eta = max(df["leading_eta"].abs().max(), df["trailing_eta"].abs().max()) 
print(f"Maximum eta recorded: {max_eta}")

print()
