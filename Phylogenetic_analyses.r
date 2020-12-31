library(diagram)
library(lubridate)
library(seraphim)
library(treeio)

writingFiles = FALSE
showingPlots = FALSE

wd1 = "A_on_IQTREE-TreeTime"
wd2 = "B_integrated_analyses"

analysis = "TreeTime_26102020"

metadata_NS_sequences = read.csv("NS_sequences_data.csv", head=T)
metadata_NY_sequences = read.csv("NY_sequences_data.csv", head=T)

NY_extent = extent(-74.27, -71.81, 40.49, 41.20)
country = crop(shapefile("NY_state_all_shapefiles/GADM_USA_0.shp"), NY_extent)
NYstate = shapefile("NY_state_all_shapefiles/State_shoreline.shp")
NYcounties = shapefile("NY_state_all_shapefiles/County_polygons.shp")
NYboroughs = shapefile("NY_state_all_shapefiles/Boroughs_NYC.shp")
NYstate = crop(spTransform(NYstate, country@proj4string), NY_extent)
NYcounties = crop(spTransform(NYcounties, country@proj4string), NY_extent)
col_boroughs = c("#60A54C75","#C3874D75","#DC442875","#4776BA75","#F9A42175")

# A. ANALYSES BASED ON A IQTREE-TREETIME TREE

	# A1. Gathering all NYC sequences in a unique fasta file

data = read.csv(paste0(wd1,"/Provided_metadata.csv"), sep=",")
dates = dmy(gsub("\\/","-",as.character(data[,"Clinical.Collection.Date"])))
files = list.files(paste0(wd1,"/All_original_sequences")); buffer = c()
for (i in 1:length(files))
	{
		buffer = c(buffer, scan(paste0(wd1,"/All_original_sequences/",files[i]), what="", sep="\n", quiet=T, blank.lines.skip=F))
	}
if (writingFiles) write(buffer, paste0(wd1,"/All_NY_sequences.fasta"))

	# A2. Preparing the input files for the discrete phylogeographic analyses 

analysisWithBoroughs = FALSE; analysisWithBoroughs = TRUE
tree = read.tree(paste0(wd1,"/",analysis,".tre")); seqIDs = tree$tip.label
locations = rep(NA, length(seqIDs)); collectionDates = rep(NA, length(seqIDs))
for (i in 1:length(seqIDs))
	{
		if (grepl("hCoV-19",seqIDs[i]))
			{
				locations[i] = unlist(strsplit(seqIDs[i],"\\/"))[2]
			}	else	{
				locations[i] = unlist(strsplit(seqIDs[i],"\\/"))[1]
			}
		if (length(unlist(strsplit(seqIDs[i],"\\|"))) == 3)
			{
				collectionDates[i] = unlist(strsplit(seqIDs[i],"\\|"))[length(unlist(strsplit(seqIDs[i],"\\|")))]
			}
		if (length(unlist(strsplit(seqIDs[i],"\\|"))) == 4)
			{
				collectionDates[i] = unlist(strsplit(seqIDs[i],"\\|"))[length(unlist(strsplit(seqIDs[i],"\\|")))-1]
			}
	}
for (i in 1:length(seqIDs))
	{
		if (seqIDs[i]%in%metadata_NY_sequences[,"GISAID.Virus.Name"])
			{
				index = which(metadata_NY_sequences[,"GISAID.Virus.Name"]==seqIDs[i])
				if (metadata_NY_sequences[index,"Zip.State"]=="NY")
					{
						locations[i] = "NY"
						if (analysisWithBoroughs == TRUE)
							{
								if (metadata_NY_sequences[index,"Zip.County"]=="New York") locations[i] = "NY-Manhattan"
								if (metadata_NY_sequences[index,"Zip.County"]=="Kings") locations[i] = "NY-Brooklyn"
								if (metadata_NY_sequences[index,"Zip.County"]=="Queens") locations[i] = "NY-Queens"
								if (metadata_NY_sequences[index,"Zip.County"]=="Bronx") locations[i] = "NY-Bronx"
								if (metadata_NY_sequences[index,"Zip.County"]=="Richmond") locations[i] = "NY-StatenIsland"
							}
					}	else	{
						locations[i] = "USA"
					}
				date = unlist(strsplit(metadata_NY_sequences[index,"Clinical.Collection.Date"],"\\/"))
				collectionDates[i] = paste(date[3],date[2],date[1],sep="-")
			}
	}
tab = cbind(seqIDs,locations,collectionDates); colnames(tab) = c("Strain","Location","Collection Data"); txt = c()
if (writingFiles) write.csv(tab, paste0(wd1,"/",analysis,".csv"), row.names=F, quote=F)
data = read.csv(paste0(wd1,"/",analysis,".csv")); tab = data[,1:3]; colnames(tab) = c("trait","location","collection_date")
for (i in 1:dim(tab)[1])
	{
		if (!grepl("NY",tab[i,"location"])) tab[i,"location"] = "other"
		txt = c(txt, paste0(">",tab[i,"trait"]),"NNNN")
	}
if (writingFiles)
	{
		if (analysisWithBoroughs != TRUE) write.table(tab, paste0(wd1,"/",analysis,".txt"), row.names=F, quote=F, sep="\t")
		if (analysisWithBoroughs == TRUE) write.table(tab, paste0(wd1,"/BSSVS_&_tips_swapping/",analysis,".txt"), row.names=F, quote=F, sep="\t")
	}
if (writingFiles) write(txt, paste0(wd1,"/",gsub("2020","20",analysis),".fasta"))

	# A3. Analysing the outputs of the preliminary discrete phylogeographic analysis 

burnIn = 1001; computingHPDInterval = FALSE # N.B.: long analysis
if (computingHPDInterval)
	{
		trees = scan(paste0(wd1,"/",gsub("2020","20",analysis),".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		indices1 = which(!grepl("tree STATE_",trees)); indices2 = which(grepl("tree STATE_",trees))
		NYstateBranches_list = rep(NA,length(trees))
		NYstateIntroductions_list = rep(NA,length(trees))
		NYstateTipBranches_list = rep(NA,length(trees))
		for (i in (burnIn+1):length(indices2))
			{
				tree1 = trees[c(indices1[1:(length(indices1)-1)],indices2[i],indices1[length(indices1)])]
				write(tree1, paste0(wd1,"/",gsub("2020","20",analysis),"_sampled_tree_",i,".tree"))
				tree2 = readAnnotatedNexus(paste0(gsub("2020","20",analysis),"_sampled_tree_",i,".tree"))
				NYstateBranches = 0; NYstateIntroductions = 0; NYstateTipBranches = 0
				for (j in 1:dim(tree2$edge)[1])
					{
						if (grepl("NY",tree2$annotations[[j]]$location))
							{
								NYstateBranches = NYstateBranches + 1
								index = which(tree2$edge[,2]==tree2$edge[j,1])
								if (grepl("NY",tree2$annotations[[index]]$location))
									{
										NYstateIntroductions = NYstateIntroductions + 1
									}
								if (!tree2$edge[j,2]%in%tree2$edge[,1])
									{
										NYstateTipBranches = NYstateTipBranches + 1
									}
							}
					}
				NYstateBranches_list[i] = NYstateBranches
				NYstateIntroductions_list[i] = NYstateIntroductions
				NYstateTipBranches_list[i] = NYstateTipBranches
				file.remove(paste0(wd1,"/",gsub("2020","20",analysis),"_sampled_tree_",i,".tree"))
			}
		quantiles = quantile(NYstateIntroductions_list[!is.na(NYstateIntroductions_list)],probs=c(0.025,0.975))
		cat("A minimum number of ",median(NYstateIntroductions_list[!is.na(NYstateIntroductions_list)])," lineage introductions (95% HPD interval = [",
			quantiles[1],"-",quantiles[2],"])"," identified from the global phylogenetic analysis of ",NYstateTipBranches," SARS-CoV-2 sampled in NY",sep="")
		# A minimum number of 116 lineage introductions (95% HPD interval = [107-127]) identified from the global phylogenetic analysis of 828 SARS-CoV-2 sampled in NY
	}

tree = readAnnotatedNexus(paste0(wd1,"/",gsub("2020","20",analysis),".tree"))
if (showingPlots)
	{
		tab = read.csv(paste0(wd1,"/",analysis,".csv"), head=T)
		samplingDates = decimal_date(ymd(gsub("\\/","-",tab[,"Collection.Data"]))); mostRecentSamplingYear = max(samplingDates, na.rm=T)
		selectedDates = decimal_date(ymd(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01")))
		rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingYear-rootHeight
		selectedLabels = c("01-01-2020","01-02-2020","01-03-2020","01-04-2020","01-05-2020")
		cols = rep("gray30",dim(tree$edge)[1]); lwds = rep(0.1,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="NY") & (tree$annotations[[i]]$location=="NY"))
							{
								cols[i] = "#57068C"; lwds[i] = 0.4
							}
					}
			}
		pdf("Figure_1_NEW.pdf", width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		# dev.new(width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("NY",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.3, col="#57068C")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.3, col="gray30", lwd=0.5)
					}
				if (tree$annotations[[i]]$location == "NY")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if (tree$annotations[[index]]$location != "NY")
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="#57068C")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
		cols = rep("gray50",dim(tree$edge)[1]); lwds = rep(0.05,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="NY") & (tree$annotations[[i]]$location=="NY"))
							{
								cols[i] = "#57068C"; lwds[i] = 0.4
							}
					}
			}
		dev.off()
	}
NYstateBranches = c(); NYstateIntroductions = c()
NYstateTipBranches = c(); sampledSequences = c()
for (i in 1:dim(tree$edge)[1])
	{
		if (tree$annotations[[i]]$location == "NY")
			{
				NYstateBranches = c(NYstateBranches,i)
				index = which(tree$edge[,2]==tree$edge[i,1])
				if (tree$annotations[[index]]$location != "NY")
					{
						NYstateIntroductions = c(NYstateIntroductions, i)
					}
				if (!tree$edge[i,2]%in%tree$edge[,1])
					{
						NYstateTipBranches = c(NYstateTipBranches, i)
						sampledSequences = c(sampledSequences, tree$tip.label[tree$edge[i,2]])
					}
			}
	}
for (i in 1:length(NYstateIntroductions))
	{
		if (i == 1) clusters1 = list()
		if (tree$edge[NYstateIntroductions[i],2]%in%tree$edge[,1])
			{
				subtree = tree_subset(tree, tree$edge[NYstateIntroductions[i],2], levels_back=0)
				clusters1[[i]] = gsub("'","",subtree$tip.label)
			}	else		{
				clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[NYstateIntroductions[i],2]])
			}
	}
for (i in 2:length(clusters1))
	{
		for (j in 1:(i-1))
			{
				if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]]))
					{
						clusters1[[j]] = clusters1[[j]][which(!clusters1[[j]]%in%clusters1[[i]])]
					}
				if (sum(clusters1[[j]]%in%clusters1[[i]]) == length(clusters1[[j]]))
					{
						clusters1[[i]] = clusters1[[i]][which(!clusters1[[i]]%in%clusters1[[j]])]
					}
			}
	}
sampledSequences = gsub("'","",sampledSequences)
zipCodes = shapefile("NY_state_all_shapefiles/ZipCodes_NY.shp")
zipCodes = shapefile("NY_state_all_shapefiles/ZipCodes_US.shp")
if (!file.exists(paste0(wd1,"/Sampling_NY_state.csv")))
	{
		data = read.csv(paste0(wd1,"/NY_sequences_data.csv"), head=T)
		samplingData = matrix(nrow=length(sampledSequences), ncol=5)
		colnames(samplingData) = c("sequenceID","collectionDate","zipCode","longitude","latitude")
		samplingData[,"sequenceID"] = sampledSequences
		for (i in 1:dim(samplingData)[1])
			{
				index = which(data[,"GISAID.Virus.Name"]==samplingData[i,"sequenceID"])
				date = dmy(gsub("\\/","-",data[index,"Clinical.Collection.Date"]))
				samplingData[i,"collectionDate"] = decimal_date(date)
				ID = unlist(strsplit(samplingData[i,"sequenceID"],"\\/"))[3]
				samplingData[i,"zipCode"] = data[index,"Zip"]
				indices = which(zipCodes@data[,"ZCTA5CE10"]==data[index,"Zip"])
				if (length(indices) > 0)
					{
						maxArea = 0; polIndex1 = 0; polIndex2 = 0
						for (j in 1:length(indices))
							{
								for (k in 1:length(zipCodes@polygons[[indices[j]]]@Polygons))
									{
										if (maxArea < zipCodes@polygons[[indices[j]]]@Polygons[[k]]@area)
											{
												maxArea = zipCodes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
											}
									}
							}
						pol = zipCodes@polygons[[polIndex1]]@Polygons[[polIndex2]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol = sps; proj4string(pol) = zipCodes@proj4string
						samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
						samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
					}
			}
		write.csv(samplingData, paste0(wd1,"/Sampling_NY_state.csv"), quote=F, row.names=F)
	}	
samplingData = read.csv(paste0(wd1,"/Sampling_NY_state.csv"), head=T)
for (i in 1:length(NYstateIntroductions))
	{
		tab = c()
		if (i == 1)
			{
				clusters2 = list(); centroids = list()
			}
		for (j in 1:length(clusters1[[i]]))
			{
				index = which(samplingData[,"sequenceID"]==clusters1[[i]][j])
				if (length(index) == 1)
					{
						line = cbind(as.numeric(samplingData[index,"collectionDate"]),as.numeric(samplingData[index,"longitude"]),as.numeric(samplingData[index,"latitude"]))
						row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
					}
			}
		colnames(tab) = c("collectionDate","longitude","latitude"); clusters2[[i]] = tab
		centroids[[i]] = cbind(mean(tab[!is.na(tab[,"longitude"]),"longitude"]), mean(tab[!is.na(tab[,"latitude"]),"latitude"]))
	}
clusterSizes = rep(NA, length(clusters1))
collectionDates = c()
for (i in 1:length(clusters1))
	{
		clusterSizes[i] = dim(clusters2[[i]])[1]
		collectionDates = c(collectionDates, clusters2[[i]][,"collectionDate"])
	}
if (showingPlots)
	{
		collectionDates_filetered = collectionDates
		dev.new(width=3.3, height=8); par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(2,2,1,1), lwd=0.2, col="gray30")
		hist(clusterSizes, breaks=50, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		hist(collectionDates_filetered, breaks=65, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=decimal_date(ymd(c("2020-02-01","2020-03-01","2020-04-01","2020-05-01"))),
			 labels=c("01-02-2020","01-03-2020","01-04-2020","01-05-2020"))
	}

	# A4. Subsampling New York broroughs with a phylogenetic clustering approach (not used)

if (!file.exists(paste0(wd1,"/BSSVS_&_tips_swapping/","TreeTime_with1000.trees")))
	{
		allTrees = scan(paste0(wd1,"/BSSVS_&_tips_swapping/",gsub("2020","20",analysis),".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		tractionFiles = 1000; index1 = which(allTrees=="\t\t;")[2]
		index2 = which(grepl("tree STATE_500000 ",allTrees))
		index3 = which(grepl("tree STATE_5000000 ",allTrees))
		interval = (index3-index2)/nberOfExtractionFiles
		selectedTrees = allTrees[c(1:index1,seq(index2+interval,index3,interval))]
		write(c(selectedTrees,"End;"), paste0(wd1,"/BSSVS_&_tips_swapping/","TreeTime_with1000.trees"))
	}
identifyingPhylogeneticClusters = FALSE
if (identifyingPhylogeneticClusters = TRUE)
	{
		selectedLocations = c("NY","NY.Bronx","NY.Brooklyn","NY.Manhattan","NY.Queens","NY.StatenIsland")
		tree = readAnnotatedNexus(paste0(wd1,"/BSSVS_&_tips_swapping/TreeTime_with1000.tree"))
		metadata = read.table(paste0(wd1,"/BSSVS_&_tips_swapping/",analysis,".txt"), head=T)
		subTrees = subtrees(tree, wait=F); sequencesToRemove = c()
		c3 = 0; clusters3 = list(); c4 = 0; clusters4 = list()
		for (i in 2:length(subTrees)) # the 1st subtree is the entire one
			{
				subTree = subTrees[i][[1]]; labels = subTree$tip.label; locations = c()
				for (j in 1:length(labels)) locations = c(locations, metadata[which(metadata[,"trait"]==labels[j]),"location"])
				if ((length(unique(locations)) == 1)&&(unique(locations)%in%selectedLocations)) { c3 = c3+1; clusters3[[c3]] = labels }
			}
		for (i in 1:length(clusters3))
			{
				nested = FALSE
				for (j in 1:length(clusters3))
					{
						if (i != j)
							{
								allSequencesIncluded = TRUE
								for (k in 1:length(clusters3[[i]]))
									{
										if (!clusters3[[i]][k]%in%clusters3[[j]]) allSequencesIncluded = FALSE
									}
								if (allSequencesIncluded == TRUE) nested = TRUE
							}
					}
				if (nested == FALSE) { c4 = c4 + 1; clusters4[[c4]] = clusters3[[i]] }
			}
		for (i in 1:length(clusters4)) sequencesToRemove = c(sequencesToRemove, sample(clusters4[[i]],length(clusters4[[i]])-1,replace=F))
	}

	# A5. Analysing the output of the BSSVS and corresponding tip swapping analyses

		# - in BEAUti, in the "States" panel: "reconstruct state change counts" and "reconstruct complete change history on tree"
		# - in the XML, add a "completeJumpHistory" to log the Markov jumps
		# - after thr run, the Markov jumps can be extracted with the Pearl script "Getting_Markov_jumps.pl":
		#	> perl Getting_Markov_jumps.pl 0 < TreeTime_261020_3.log > TreeTime_261020_3.txt

locations = c("NY","NY-Bronx","NY-Brooklyn","NY-Manhattan","NY-Queens","NY-StatenIsland","other")
if (!file.exists(paste0(wd1,"/BSSVS_&_tips_swapping/TreeTime_261020_3.txt")))
	{
		wd = getwd(); setwd(paste0(wd,"/",wd1,"/BSSVS_&_tips_swapping/"))
		system("perl Getting_Markov_jumps.pl 0 < TreeTime_261020_3.log > TreeTime_261020_3.txt"); setwd(wd)
	}
if (!file.exists("BSSVS_MarkovJumps1.csv"))
	{
		wd = getwd(); setwd(paste0(wd,"/",wd1,"/BSSVS_&_tips_swapping/"))
		tab = read.table("TreeTime_261020_3.txt", head=T)
		burnIn = 1000000; tab = tab[which(tab[,"state"]>burnIn),]
		MJs = matrix(nrow=length(locations), ncol=length(locations))
		states = unique(tab[,"state"])
		for (i in 1:length(locations))
			{
				for (j in 1:length(locations))
					{
						sub = tab[which((tab[,"from"]==locations[i])&(tab[,"to"]==locations[j])),]
						MJs_post = c(); # print(c(i,j))
						for (k in 1:length(states))
							{
								MJ = which(sub[,"state"]==states[k])
								MJs_post = c(MJs_post, length(MJ))
							}
						MJs[i,j] = mean(MJs_post) # median ??
					}
			}
		row.names(MJs) = locations; colnames(MJs) = locations
		write.table(round(MJs,1), "BSSVS_MarkovJumps1.csv", quote=F, sep=",")
		setwd(wd)
	}	else		{
		MJs = read.csv(paste0(wd1,"/BSSVS_&_tips_swapping/BSSVS_MarkovJumps1.csv"), header=T)
	}
locations = c("NY","NY.Bronx","NY.Brooklyn","NY.Manhattan","NY.Queens","NY.StatenIsland","other")
log1 = read.table(paste0(wd1,"/BSSVS_&_tips_swapping/",gsub("2020","20",analysis),"_1.log"), header=T)
log2 = read.table(paste0(wd1,"/BSSVS_&_tips_swapping/",gsub("TreeTime","TipsSwap",gsub("2020","20",analysis)),"_1.log"), header=T)
log1 = log1[(round(dim(log1)[1]/10)+1):dim(log1)[1],]
log2 = log2[(round(dim(log2)[1]/10)+1):dim(log2)[1],]
rates = matrix(nrow=length(locations), ncol=length(locations))
BFs1 = matrix(nrow=length(locations), ncol=length(locations))
BFs2 = matrix(nrow=length(locations), ncol=length(locations))
row.names(rates) = locations; colnames(rates) = locations
row.names(BFs1) = locations; colnames(BFs1) = locations
row.names(BFs2) = locations; colnames(BFs2) = locations
for (i in 1:length(locations))
	{
		for (j in 1:length(locations))
			{
				if (i != j)
					{
						colName = paste0("location.indicators.",locations[i],".",locations[j])
						index1 = which(colnames(log1)==colName); index2 = which(colnames(log2)==colName)
						p = sum(log1[,index1]==1)/dim(log1)[1]
						K = 42 # length(locations)*(length(locations)-1) # K shoulf be divided by 2 if "symetric" case
						q = (log(2)+K-1)/(K*(K-1))
						BFs1[i,j] = (p/(1-p))/(q/(1-q))
						p1 = sum(log1[,index1]==1)/dim(log1)[1]
						p2 = sum(log2[,index2]==1)/dim(log2)[1]
						BFs2[i,j] = (p1/(1-p1))/(p2/(1-p2))
						index3 = which(colnames(log1)==paste0("location.rates.",locations[i],".",locations[j]))
						rates[i,j] = median(log1[,index3]*log1[,index1])
					}
			}
	}
MJs1 = MJs; MJs2 = MJs; rates[BFs2<3] = NA; MJs1[BFs1<3] = NA; MJs2[BFs2<3] = NA; BFs1 = round(BFs1, 1); BFs2 = round(BFs2, 1)
if (writingFiles)
	{
		write.table(BFs1, paste0(wd1,"/BSSVS_&_tips_swapping/BSSVS_BayesFactor1.csv"), sep=",", quote=F)
		write.table(BFs2, paste0(wd1,"/BSSVS_&_tips_swapping/BSSVS_BayesFactor2.csv"), sep=",", quote=F)
		write.table(rates, paste0(wd1,"/BSSVS_&_tips_swapping/BSSVS_transitionRate.csv"), sep=",", quote=F)
		write.table(MJs2, paste0(wd1,"/BSSVS_&_tips_swapping/BSSVS_MarkovJumps2.csv"), sep=",", quote=F)
	}
if (showingPlots)
	{
		pdf(paste0("Figure_B_NEW.pdf"), width=7, height=7) # dev.new(width=7, height=7)
		par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
		centroids = coordinates(NYboroughs); MJs = MJs1; MJs = MJs2
		plot(NYboroughs, col=col_boroughs, border="gray60", lwd=0.5)
		points(centroids, cex=1, pch=16, col="gray30")
		minMJ = min(MJs[MJs!=0], na.rm=T); maxMJ = max(MJs[MJs!=0], na.rm=T)
		for (i in 1:length(locations))
			{
				index1 = which(gsub(" ","",NYboroughs@data[,"boro_name"])==gsub("NY.","",locations[i]))
				for (j in 1:length(locations))
					{
						index2 = which(gsub(" ","",NYboroughs@data[,"boro_name"])==gsub("NY.","",locations[j]))
						if ((length(index1)>0)&(length(index2)>0)&(!is.na(MJs[i,j]))&(MJs[i,j]!=0))
							{
								LWD = (((MJs[i,j]-minMJ)/(maxMJ-minMJ))*8.5)+0.5; arrow = (0.5*(MJs[i,j]/maxMJ))+0.15
								curvedarrow(centroids[index1,], centroids[index2,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=-0.15, dr=NA, endhead=F, arr.type="triangle")
							}
					}
			}
		MJ = 3; LWD = (((MJ-minMJ)/(maxMJ-minMJ))*8.5)+0.5; arrow = (0.5*(MJ/maxMJ))+0.15
		curvedarrow(cbind(-74.23,40.900), cbind(-74.18,40.900), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
					lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
		MJ = 20; LWD = (((MJ-minMJ)/(maxMJ-minMJ))*8.5)+0.5; arrow = (0.5*(MJ/maxMJ))+0.15
		curvedarrow(cbind(-74.23,40.887), cbind(-74.18,40.887), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
					lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
		MJ = 50; LWD = (((MJ-minMJ)/(maxMJ-minMJ))*8.5)+0.5; arrow = (0.5*(MJ/maxMJ))+0.15
		curvedarrow(cbind(-74.23,40.874), cbind(-74.18,40.874), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
					lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
		mtext("3 Markov jumps", at=-74.17, line=-2.8, adj=0, col="gray30", cex=0.7)
		mtext("20 Markov jumps", at=-74.17, line=-3.75, adj=0, col="gray30", cex=0.7)
		mtext("50 Markov jumps", at=-74.17, line=-4.7, adj=0, col="gray30", cex=0.7)
		legend(-74.238, 40.85, NYboroughs@data[c(3,1,5,4,2),"boro_name"], pch=16, col=gsub("75","95",col_boroughs[c(3,1,5,4,2)]),
			   pt.cex=1.2, cex=0.75, y.intersp=1.2, x.intersp=1.2, text.col="gray30", bty="n")
		dev.off()
	}

	# A6. Preparing the continuous phylogeographic analyses (RRW, Cauchy model)

template = scan(paste0(wd1,"/Template_file_RRW2.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
phylogeographicRuns = paste0(wd1,"/Phylogeographic_runs"); xml = c()
if (!file.exists(paste0(phylogeographicRuns,"/All_clades_NEW.xml")))
	{
		sink(file=paste0(phylogeographicRuns,"/All_clades_NEW.xml"))
		for (i in 1:length(template))
			{
				cat(template[i],"\n")
				if (grepl("Insert taxa blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
										for (k in 1:dim(clusters2[[j]])[1])
											{
												if (!is.na(clusters2[[j]][k,"longitude"]))
													{
														cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[j]])[k],"\">","\n"))
														cat(paste0("\t\t\t<date value=\"",clusters2[[j]][k,"collectionDate"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
														cat("\t\t\t<attr name=\"latitude\">\n")
														cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"],"\n"))
														cat("\t\t\t</attr>\n")
														cat("\t\t\t<attr name=\"longitude\">\n")
														cat(paste0("\t\t\t\t",clusters2[[j]][k,"longitude"],"\n"))
														cat("\t\t\t</attr>\n")
														cat("\t\t\t<attr name=\"coordinates\">\n")
														cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"]," ",clusters2[[j]][k,"longitude"],"\n"))
														cat("\t\t\t</attr>\n")
														cat("\t\t</taxon>\n")
													}
											}
										cat("\t</taxa>","\n")
									}
							}
					}
				if (grepl("Insert alignment blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
										for (k in 1:dim(clusters2[[j]])[1])
											{
												if (!is.na(clusters2[[j]][k,"longitude"]))
													{
														cat("\t\t<sequence>\n")
														cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[j]])[k],"\"/>","\n"))
														cat("\t\t\tNNNN\n")
														cat("\t\t</sequence>\n")
													}
											}
										cat("\t</alignment>","\n")
									}
							}
					}
				if (grepl("Insert pattern blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
										cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
										cat("\t</patterns>","\n")
									}
							}
					}
				if (grepl("Insert starting tree blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										tre = tree_subset(tree, tree$edge[NYstateIntroductions[j],2], levels_back=0)
										tips = row.names(clusters2[[j]]); tips = tips[which(!is.na(clusters2[[j]][,"longitude"]))]
										tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
										if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
										write.tree(tre, paste0(phylogeographicRuns,"/Clade_",j,".tre"))
										tre = scan(paste0(phylogeographicRuns,"/Clade_",j,".tre"), what="", sep="\n", quiet=T)
										txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
										write(txt, paste0(phylogeographicRuns,"/Clade_",j,".tre"))
										cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Clade_",j,".tre\">","\n"))
										cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
										cat("\t</empiricalTreeDistributionModel>","\n")
									}
							}
					}
				if (grepl("Insert tree model blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<treeModel id=\"treeModel_",j,"\">","\n"))
										cat(paste0("\t\t<coalescentTree idref=\"startingTree_",j,"\"/>","\n"))
										cat("\t\t<rootHeight>","\n")
										cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_",j,"\"/>","\n"))
										cat("\t\t</rootHeight>","\n")
										cat("\t\t<nodeHeights internalNodes=\"true\">","\n")
										cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_",j,"\"/>","\n"))
										cat("\t\t</nodeHeights>","\n")
										cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">","\n")
										cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_",j,"\"/>","\n"))
										cat("\t\t</nodeHeights>","\n")
										cat("\t</treeModel>","\n")
									}
							}
					}
				if (grepl("Insert arbitraryBranchRates blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<arbitraryBranchRates id=\"coordinates.diffusion.branchRates",j,"\">","\n"))
										cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat("\t\t<rates>","\n")
										cat(paste0("\t\t\t<parameter id=\"coordinates.diffusion.rates",j,"\" lower=\"0.0\"/>","\n"))
										cat("\t\t</rates>","\n")
										cat("\t</arbitraryBranchRates>","\n")
									}
							}
					}
				if (grepl("Insert distributionLikelihood blocks 1",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<distributionLikelihood id=\"coordinates.diffusion.prior",j,"\">","\n"))
										cat("\t\t<data>","\n")
										cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
										cat("\t\t</data>","\n")
										cat("\t\t<distribution>","\n")
										cat(paste0("\t\t\t<onePGammaDistributionModel>","\n"))
										cat("\t\t\t\t<shape>","\n")
										cat("\t\t\t\t\t<parameter value=\"0.5\"/>","\n")
										cat("\t\t\t\t</shape>","\n")
										cat("\t\t\t</onePGammaDistributionModel>","\n")
										cat("\t\t</distribution>","\n")
										cat("\t</distributionLikelihood>","\n")
									}
							}
					}
				if (grepl("Insert coordinates.traitLikelihood blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<multivariateTraitLikelihood id=\"coordinates.traitLikelihood",j,"\" traitName=\"coordinates\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">","\n"))
										cat("\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
										cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>"))
										cat("\t\t<traitParameter>","\n")
										cat(paste0("\t\t\t<parameter id=\"leaf.coordinates",j,"\"/>","\n"))
										cat("\t\t</traitParameter>","\n")
										cat("\t\t<conjugateRootPrior>","\n")
										cat("\t\t\t<meanParameter>","\n")
										cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>","\n")
										cat("\t\t\t</meanParameter>","\n")
										cat("\t\t\t<priorSampleSize>","\n")
										cat("\t\t\t\t<parameter value=\"0.000001\"/>","\n")
										cat("\t\t\t</priorSampleSize>","\n")
										cat("\t\t</conjugateRootPrior>","\n")
										cat(paste0("\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
										cat("\t</multivariateTraitLikelihood>","\n")
									}
							}
					}
				if (grepl("Insert continuousDiffusionStatistic blocks 1",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<continuousDiffusionStatistic id=\"coordinates.diffusionRate",j,"\" greatCircleDistance=\"true\">","\n"))
										cat(paste0("\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
										cat("\t</continuousDiffusionStatistic>","\n")
									}
							}
					}
				if (grepl("Insert scaleOperator blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">","\n"))
										cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
										cat("\t\t</scaleOperator>","\n")
									}
							}
					}
				if (grepl("Insert precisionGibbsOperator blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t<precisionGibbsOperator weight=\"2\">","\n"))
										cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
										cat("\t\t\t<multivariateWishartPrior idref=\"coordinates.precisionPrior\"/>","\n")
										cat("\t\t</precisionGibbsOperator>","\n")
									}
							}
					}
				if (grepl("Insert distributionLikelihood blocks 2",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<distributionLikelihood idref=\"coordinates.diffusion.prior",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert multivariateTraitLikelihood blocks 1",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert continuousDiffusionStatistic blocks 2",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"coordinates.diffusionRate",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert multivariateTraitLikelihood blocks 2",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("<!-- Insert logTree blocks -->",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t<logTree id=\"treeFileLog",j,"\" logEvery=\"100000\" nexusFormat=\"true\" fileName=\"Clade_",j,".trees\" sortTranslationTable=\"true\">","\n"))
										cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat("\t\t\t<joint idref=\"joint\"/>","\n")
										cat("\t\t\t<trait name=\"coordinates\" tag=\"coordinates\">","\n")
										cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
										cat("\t\t\t</trait>","\n")
										cat("\t\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
										cat("\t\t\t<trait name=\"rate\" tag=\"coordinates.rate\">","\n")
										cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
										cat("\t\t\t</trait>","\n")
										cat("\t\t</logTree>","\n")
									}
							}
					}
			}
		sink(NULL)
	}

	# A7. Running BEAST and building the maximum clade consensus (MCC) tree

source("MCC_tree_extraction.r"); analyses = c()
runningNewAnalyses = FALSE; wd = getwd()
setwd(paste0(wd,"/",wd1,"/Phylogeographic_runs/"))
for (i in 1:length(clusters2))
	{
		if ((dim(clusters2[[i]])[1] >= 3)&(sum(!is.na(clusters2[[i]][,"longitude"])) >= 3)) analyses = c(analyses, paste0("Clade_",i))
	}
if (runningNewAnalyses)
	{
		system("java -jar beast_1104.jar All_clades.xml", ignore.stdout=T, ignore.stderr=F)
		for (i in 1:length(analyses))
			{
system(paste0("BEAST_1104/bin/treeannotator -burninTrees 501 -heights keep ",analyses[i],".trees ",analyses[i],".tree"), ignore.stdout=F, ignore.stderr=F)
			}
	}
setwd(wd)

	# A8. Extracting spatio-temporal information embedded in MCC and posterior trees

setwd(paste0(wd,"/",wd1,"/Phylogeographic_runs/"))
for (i in 1:length(analyses))
	{
		if (!file.exists(paste0(analyses[i],".csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
				mcc_tre = readAnnotatedNexus(paste0(analyses[i],".tree")); dates = c()
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
				mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
				write.csv(mcc_tab, paste0(analyses[i],".csv"), row.names=F, quote=F)
			}
	}
nberOfTreesToSample = 1000; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 5
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0(analyses[i],"_ext")
		if (!file.exists(paste0(localTreesDirectory,"/TreeExtractions_1.csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2]); burnIn = 501
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
				allTrees = scan(file=paste0(analyses[i],".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
			}
	}
cladesToExclude = c()
if (!file.exists("All_clades.csv"))
	{
		for (i in 1:length(analyses))
			{
				tab = read.csv(paste0(analyses[i],".csv"), head=T)
				if (i == 1)
					{
						all = tab
					}	else	{
						if (!analyses[i]%in%cladesToExclude)
							{
								maxNodeID = max(all[,c("node1","node2")])
								tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
								all = rbind(all, tab)
							}
					}
			}
		write.csv(all, "All_clades.csv", row.names=F, quote=F)
	}
dir.create(file.path("All_clades_ext"), showWarnings=F)
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:nberOfExtractionFiles)
	{
		if (!file.exists(paste0("All_clades_ext/TreeExtractions_",i,".csv")))
			{
				for (j in 1:length(analyses))
					{
						tab = read.csv(paste0(analyses[j],"_ext/TreeExtractions_",i,".csv"), head=T)
						if (j == 1)
							{
								all = tab
							}	else	{
								if (!analyses[i]%in%cladesToExclude)
									{
										maxNodeID = max(all[,c("node1","node2")])
										tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
										all = rbind(all, tab)
									}
							}
					}
				write.csv(all, paste0("All_clades_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
	}
mcc = read.csv("All_clades.csv")
tmp = mcc[which(mcc[,"node2"]%in%mcc[,"node1"]),]
write.csv(tmp, "Internals.csv", row.names=F, quote=F)
dir.create(file.path("All_clades_int"), showWarnings=F)
for (i in 1:nberOfExtractionFiles)
	{
		if (!file.exists(paste0("All_clades_int/TreeExtractions_",i,".csv")))
			{
				tab = read.csv(paste0("All_clades_ext/TreeExtractions_",i,".csv"), head=T)
				tmp = tab[which(tab[,"node2"]%in%tab[,"node1"]),]
				write.csv(tmp, paste0("All_clades_int/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
	}
mcc = read.csv("All_clades.csv"); indices = c()
if (!file.exists("5_boroughs.csv"))
	{
		for (i in 1:dim(mcc)[1])
			{
				keepNode1 = FALSE; keepNode2 = FALSE
				for (j in 1:length(NYboroughs@polygons))
					{
						for (k in 1:length(NYboroughs@polygons[[j]]@Polygons))
							{
								pol = NYboroughs@polygons[[j]]@Polygons[[k]]
								if (point.in.polygon(mcc[i,"startLon"],mcc[i,"startLat"],pol@coords[,1],pol@coords[,2])) keepNode1 = TRUE
								if (point.in.polygon(mcc[i,"endLon"],mcc[i,"endLat"],pol@coords[,1],pol@coords[,2])) keepNode2 = TRUE
							}
					}
				if ((keepNode1 == TRUE)&(keepNode2 == TRUE)) indices = c(indices, i)
			}
		write.csv(mcc[indices,], "5_boroughs.csv", row.names=F, quote=F)
	}
dir.create(file.path("5_boroughs"), showWarnings=F)
for (i in 1:nberOfExtractionFiles)
	{
		if (!file.exists(paste0("5_boroughs/TreeExtractions_",i,".csv")))
			{
				tab = read.csv(paste0("All_clades_ext/TreeExtractions_",i,".csv"), head=T); indices = c()
				for (j in 1:dim(tab)[1])
					{
						keepNode1 = FALSE; keepNode2 = FALSE
						for (k in 1:length(NYboroughs@polygons))
							{
								for (l in 1:length(NYboroughs@polygons[[k]]@Polygons))
									{
										pol = NYboroughs@polygons[[k]]@Polygons[[l]]
										if (point.in.polygon(tab[j,"startLon"],tab[j,"startLat"],pol@coords[,1],pol@coords[,2])) keepNode1 = TRUE
										if (point.in.polygon(tab[j,"endLon"],tab[j,"endLat"],pol@coords[,1],pol@coords[,2])) keepNode2 = TRUE
									}	
							}
						if ((keepNode1 == TRUE)&(keepNode2 == TRUE)) indices = c(indices, j)
					}
				write.csv(tab[indices,], paste0("5_boroughs/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
	}
setwd(wd)

	# A9. Generating a dispersal history graph (mapped MCC trees, 80% HPD polygons)

if (showingPlots)
	{
		phylogeographicRuns = paste0(wd1,"/Phylogeographic_runs")
		localTreesDirectory1 = paste0(phylogeographicRuns,"/All_clades_ext")
		localTreesDirectory2 = paste0(phylogeographicRuns,"/5_boroughs")
		mcc1 = read.csv(paste0(phylogeographicRuns,"/All_clades.csv"), head=T); startDatum = min(mcc1[,"startYear"])
		mcc2 = read.csv(paste0(phylogeographicRuns,"/5_boroughs.csv"), head=T); startDatum = min(mcc1[,"startYear"])
		nberOfExtractionFiles = 1000; percentage = 80; prob = percentage/100; precision = 1/(365/7)
		MRCApolygons = list(); MRCAdates = list(); MRCAs = which(!mcc1[,"node1"]%in%mcc1[,"node2"])
		colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(141)[16:116])
		for (i in 1:length(MRCAs))
			{
				nodeCoordinates = matrix(nrow=nberOfExtractionFiles, ncol=2); colnames(nodeCoordinates) = c("lon","lat")
				for (j in 1:nberOfExtractionFiles)
					{
						tab = read.csv(paste0(phylogeographicRuns,"/All_clades_ext/TreeExtractions_",j,".csv"), head=T)
						nodeCoordinates[j,"lon"] = tab[MRCAs[i],"startLon"]; nodeCoordinates[j,"lat"] = tab[MRCAs[i],"startLat"]
					}
				H = Hpi(cbind(nodeCoordinates[,"lon"],nodeCoordinates[,"lat"]))
				kde = kde(cbind(nodeCoordinates[,"lon"],nodeCoordinates[,"lat"]), H=H, compute.cont=T, gridsize=c(1000,1000))
				contourLevel = contourLevels(kde, prob=(1-(percentage/100))); polygons = list()
				contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
				for (j in 1:length(contourLines)) polygons[[j]] = Polygon(cbind(contourLines[[j]]$x,contourLines[[j]]$y))
				ps = Polygons(polygons,1); contourPolygons = SpatialPolygons(list(ps))
				spdf = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
				names(spdf) = round(tab[MRCAs[i],"startYear"],3); MRCApolygons[[i]] = spdf
			}
		minYear = min(mcc1[,"startYear"]); maxYear = max(mcc1[,"endYear"])
		MRCApolygons_colours = rep(NA, length(MRCApolygons))
		for (i in 1:length(MRCApolygons))
			{
				date = as.numeric(names(MRCApolygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				MRCApolygons_colours[i] = paste0(colourScale[polygon_index],"40")
			}
		polygons1 = suppressWarnings(spreadGraphic2(localTreesDirectory1, nberOfExtractionFiles, prob, startDatum, precision))
		startYears_indices1 = (((mcc1[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices1 = (((mcc1[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours1 = colourScale[startYears_indices1]; endYears_colours1 = colourScale[endYears_indices1]
		polygons_colours1 = rep(NA, length(polygons1))
		for (i in 1:length(polygons1))
			{
				date = as.numeric(names(polygons1[[i]])); polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours1[i] = paste0(colourScale[polygon_index],"40")
			}
		polygons2 = suppressWarnings(spreadGraphic2(localTreesDirectory2, nberOfExtractionFiles, prob, startDatum, precision))
		startYears_indices2 = (((mcc2[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices2 = (((mcc2[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours2 = colourScale[startYears_indices2]; endYears_colours2 = colourScale[endYears_indices2]
		polygons_colours2 = rep(NA, length(polygons2))
		for (i in 1:length(polygons2))
			{
				date = as.numeric(names(polygons2[[i]])); polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours2[i] = paste0(colourScale[polygon_index],"40")
			}
	}
if (showingPlots)
	{
		selectedDates = c(minYear,decimal_date(ymd(c("2020-03-01","2020-04-01","2020-05-01"))),maxYear)
		selectedLabels = c("","01-03-2020","01-04-2020","01-05-2020",""); dev.new(width=4, height=2)
		par(mar=c(3,3,1,1), oma=c(0,0,0,0), lwd=0.2, col="gray30", lwd=0.2, cex=0.6)
		plot(density(mcc1[which(!mcc[,"node2"]%in%mcc[,"node1"]),"endYear"]), xlim=c(minYear,maxYear), main=NA, xlab=NA, cex.lab=0.7, axes=F, col="gray30", frame=F)
		axis(1, cex=0.6, lwd=0.2, lwd.tick=0.2, tck=-0.025, col="gray30", line=0, mgp=c(0,0.5,0), at=selectedDates, labels=selectedLabels)
		axis(2, cex=0.6, lwd=0.2, lwd.tick=0.2, tck=-0.025, col="gray30", line=0, mgp=c(0,0.5,0), ylab="density")
	}
if (showingPlots)
	{
		tab = read.csv("NY_epidemiological_data/NY_state_COVID_data.csv", head=T)
		specimen_dates = gsub("2270-","2020-",specimen_dates); specimen_dates = gsub("2271-","2020-",specimen_dates)
		specimen_dates = gsub("2272-","2020-",specimen_dates); specimen_dates = gsub("2273-","2020-",specimen_dates)
		specimen_dates = gsub("2297-","2020-",specimen_dates); specimen_dates = gsub("2299-","2020-",specimen_dates)
		specimen_dates = unique(tab[,"specimen_date"]); hospitalisationDays = c()
		for (i in 1:length(specimen_dates))
			{
				indices = which(tab[,"specimen_date"]==specimen_dates[i])
				extract_dates = decimal_date(mdy(tab[indices,"extract_date"]))
				count = tab[indices[which(extract_dates==max(extract_dates))],"Number_hospitalized"]
				hospitalisationDays = c(hospitalisationDays, rep(decimal_date(mdy(specimen_dates[i])),count))
			}
		selectedDates = c(minYear,decimal_date(ymd(c("2020-03-01","2020-04-01","2020-05-01"))),maxYear)
		selectedLabels = c("","01-03-2020","01-04-2020","01-05-2020",""); dev.new(width=4, height=2)
		par(mar=c(3,3,1,1), oma=c(0,0,0,0), lwd=0.2, col="gray30", lwd=0.2, cex=0.6)
		hist(hospitalisationDays, breaks=275, xlim=c(minYear,maxYear), main=NA, xlab=NA, cex.lab=0.7, axes=F, col="gray90")
		axis(1, cex=0.6, lwd=0.2, lwd.tick=0.2, tck=-0.025, col="gray30", line=0, mgp=c(0,0.5,0), at=selectedDates, labels=selectedLabels)
		axis(2, cex=0.6, lwd=0.2, lwd.tick=0.2, tck=-0.025, col="gray30", line=0, mgp=c(0,0.5,0), ylab="density")
	}
if (showingPlots)
	{
		figure = "Ca"; plotState = T; plotCladeOfSize2 = T; plotOnlyMRCANodes = F; plotOnlyInternalNodes = F; plotOnlyTipNodes = F; plotPolygons = T; croppingPolygons = F
		figure = "Cb"; plotState = T; plotCladeOfSize2 = F; plotOnlyMRCANodes = F; plotOnlyInternalNodes = T; plotOnlyTipNodes = F; plotPolygons = T; croppingPolygons = F
		figure = "Cc"; plotState = T; plotCladeOfSize2 = T; plotOnlyMRCANodes = F; plotOnlyInternalNodes = F; plotOnlyTipNodes = T; plotPolygons = T; croppingPolygons = F
		figure = "Da"; plotState = F; plotCladeOfSize2 = T; plotOnlyMRCANodes = F; plotOnlyInternalNodes = F; plotOnlyTipNodes = F; plotPolygons = T; croppingPolygons = F
		figure = "Db"; plotState = F; plotCladeOfSize2 = F; plotOnlyMRCANodes = F; plotOnlyInternalNodes = T; plotOnlyTipNodes = F; plotPolygons = T; croppingPolygons = F
		figure = "Dc"; plotState = F; plotCladeOfSize2 = T; plotOnlyMRCANodes = F; plotOnlyInternalNodes = F; plotOnlyTipNodes = T; plotPolygons = T; croppingPolygons = F
		cexNode = 0.7; LWD = 1.0
		if (plotState == TRUE)
			{
				pdf(paste0("Figure_",figure,"_NEW.pdf"), width=11.0, height=4.5) # dev.new(width=11.0, height=3.9)
				par(oma=c(0,0,0,0), mar=c(0,1,0,1), lwd=0.2, col="gray30")
				plot(as(NY_extent, "SpatialPolygons"), border=NA)
				plot(country, col="gray95", border=NA, add=T)		
				plot(NYstate, col="gray90", border=NA, add=T)
				plot(as(NY_extent, "SpatialPolygons"), add=T, border="gray30")
			}	else	{
				pdf(paste0("Figure_",figure,"_NEW.pdf"), width=7, height=7) # dev.new(width=7, height=7)
				par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
				plot(NYboroughs, col="gray90", border="white", lwd=0.5)
			}
		if (plotPolygons == TRUE)
			{
				if (plotState == TRUE)
					{
						polygons = polygons1; polygons_colours = polygons_colours1
					}
				if (plotState != TRUE)
					{
						polygons = polygons2; polygons_colours = polygons_colours2
					}
				if (plotOnlyMRCANodes == TRUE) polygons = MRCApolygons
				for (i in length(polygons):1)
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]
						if (croppingPolygons == TRUE)
							{
								if (plotState == TRUE) pol = crop(pol, NYstate)
								if (plotState != TRUE) pol = crop(pol, NYboroughs)
							}
						plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
					}
			}
		if (plotState == TRUE)
			{
				mcc = mcc1; startYears_colours = startYears_colours1; endYears_colours = endYears_colours1
			}
		if (plotState != TRUE)
			{
				mcc = mcc2; startYears_colours = startYears_colours2; endYears_colours = endYears_colours2
			}
		selectedBranches = 1:dim(mcc)[1]
		if (plotOnlyMRCANodes == TRUE) selectedBranches = selectedBranches[which(!mcc[,"node1"]%in%mcc[,"node2"])]
		if (plotOnlyInternalNodes == TRUE) selectedBranches = selectedBranches[which(mcc[,"node2"]%in%mcc[,"node1"])]
		if (plotOnlyTipNodes == TRUE) selectedBranches = selectedBranches[which(!mcc[,"node2"]%in%mcc[,"node1"])]
		if (plotOnlyMRCANodes != TRUE)
			{
				for (i in selectedBranches)
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
								    arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
			}
		if (plotCladeOfSize2 == TRUE)
			{
				for (i in 1:length(clusters2))
					{
						if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
							{
								if (sum(!is.na(clusters2[[i]][,"longitude"])) == 2)
									{
										indices = which(!is.na(clusters2[[i]][,"longitude"]))
										if (length(indices) == 2)
											{
												if (plotState == TRUE)
													{
														buffer = indices
													}	else	{
														keepNode1 = FALSE; keepNode2 = FALSE; buffer = c()
														for (j in 1:length(NYboroughs@polygons))
															{
																for (k in 1:length(NYboroughs@polygons[[j]]@Polygons))
																	{
																		pol = NYboroughs@polygons[[j]]@Polygons[[k]]
										if (point.in.polygon(clusters2[[i]][indices[1],"longitude"],clusters2[[i]][indices[1],"latitude"],pol@coords[,1],pol@coords[,2])) keepNode1 = TRUE
										if (point.in.polygon(clusters2[[i]][indices[2],"longitude"],clusters2[[i]][indices[2],"latitude"],pol@coords[,1],pol@coords[,2])) keepNode1 = TRUE
																	}	
															}
														if ((keepNode1 == TRUE)&(keepNode2 == TRUE)) buffer = c(buffer, j)
													}
												if (length(buffer) == 2)
													{
														curvedarrow(cbind(clusters2[[i]][buffer[1],"longitude"],clusters2[[i]][buffer[1],"latitude"]),
																	cbind(clusters2[[i]][buffer[2],"longitude"],clusters2[[i]][buffer[2],"latitude"]),
																	arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="gray10", arr.col=NA, arr.pos=F,
																	curve=0.1, dr=NA, endhead=F)
													}
											}
									}
							}
					}
			}
		if ((plotOnlyMRCANodes != TRUE)&(plotOnlyInternalNodes != TRUE))
			{
				for (i in 1:length(clusters2))
					{
						if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
							{
								for (j in 1:dim(clusters2[[i]])[1])
									{
										if (!is.na(clusters2[[i]][j,"longitude"]))
											{
												plotTheNode = FALSE
												if (plotState == TRUE)
													{
														plotTheNode = TRUE
													}	else	{
														for (k in 1:length(NYboroughs@polygons))
															{
																for (l in 1:length(NYboroughs@polygons[[k]]@Polygons))
																	{
																		pol = NYboroughs@polygons[[k]]@Polygons[[l]]
										if (point.in.polygon(clusters2[[i]][j,"longitude"],clusters2[[i]][j,"latitude"],pol@coords[,1],pol@coords[,2])) plotTheNode = TRUE
																	}	
															}
													}
												if (plotTheNode)
													{
														index = (((clusters2[[i]][j,"collectionDate"]-minYear)/(maxYear-minYear))*100)+1
														points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=16, col=colourScale[index], cex=cexNode)
														points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
													}
											}
									}
							}
					}
			}
		for (i in rev(selectedBranches))
			{
				if (plotOnlyMRCANodes != TRUE)
					{
						if (!mcc[i,"node1"]%in%mcc[selectedBranches,"node2"])
							{
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
							}
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)					
					}	else		{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=1.5*cexNode)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=1.5*cexNode, lwd=0.4)						
					}
			}
		if (plotState == TRUE)
			{
				if (plotOnlyInternalNodes != TRUE)
					{
						selectedDates = decimal_date(ymd(c("2020-03-01","2020-04-01","2020-05-01")))
						selectedLabels = c("01-03-2020","01-04-2020","01-05-2020")	
					}	else		{
						selectedDates = mcc[which(!mcc[,"node1"]%in%mcc[,"node2"]),"endYear"]
						selectedLabels = rep("", length(selectedDates))
					}
				rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.65,0.915,0.230,0.242),
					 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
				     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.00,0),
				     at=selectedDates, labels=selectedLabels))
			}
		dev.off()
	}

	# A.10. Estimating and plotting dispersal statistics associated with lineages

recomputingStatistics = FALSE
if (recomputingStatistics == TRUE)
	{
		nberOfExtractionFiles = 1000; timeSlices = 100; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 5; slidingWindow = 1/(365/14)
		localTreesDirectory = paste0(wd1,"/Phylogeographic_runs/All_clades_ext"); outputName = paste0(wd1,"/All_dispersal_statistics/",analysis)
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
		tab = read.table(paste0(wd1,"/All_dispersal_statistics/",analysis,"_estimated_dispersal_statistics.txt"), header=T)
		vS = round(tab[,"weighted_branch_dispersal_velocity"]/366,1); cat(median(vS)," km/day (95% HPD interval = [",quantile(vS,0.025),"-",quantile(vS,0.975),"])",sep="")
	}

# B. INTEGRATED ANALYSES BASED ON THE MAJOR NY CLADE

nberOfReplicates = 10

	# B1. Preparing the integrated discrete and continuous phylogeographic analyses

clusters_sizes = rep(NA, length(clusters2))
for (i in 1:length(clusters2)) clusters_sizes[i] = dim(clusters2[[i]])[1]
major_clade = clusters2[[which(clusters_sizes==max(clusters_sizes))]]
for (i in 1:dim(major_clade)[1])
	{
		index = which(metadata_NY_sequences[,"GISAID.Virus.Name"] == row.names(major_clade)[i])
		date = unlist(strsplit(metadata_NY_sequences[index,"Clinical.Collection.Date"],"\\/"))
		major_clade[i,"collectionDate"] = paste(date[3],date[2],date[1],sep="-")
	}
major_clade = cbind(row.names(major_clade), major_clade)
row.names(major_clade) = c(); colnames(major_clade)[1:2] = c("trait","collection_date")
major_clade_boroughs = matrix(nrow=dim(major_clade)[1], ncol=1)
major_clade_sequences = rep(NA, dim(major_clade)[1])
colnames(major_clade_boroughs) = "location"
for (i in 1:dim(major_clade)[1])
	{
		borough = c()
		for (j in 1:length(NYboroughs@polygons))
			{
				for (k in 1:length(NYboroughs@polygons[[j]]@Polygons))
					{
						pol = NYboroughs@polygons[[j]]@Polygons[[k]]@coords
						if (point.in.polygon(major_clade[i,"longitude"],major_clade[i,"latitude"],pol[,1],pol[,2]) == 1)
							{
								borough = c(borough, gsub(" ","",NYboroughs@data[j,"boro_name"]))
							}
					}
			}
		if (length(borough) == 1)
			{
				major_clade_boroughs[i,1] = borough
			}	else		{
				inState = FALSE
				for (j in 1:length(NYstate@polygons))
					{
						for (k in 1:length(NYstate@polygons[[j]]@Polygons))
							{
								pol = NYstate@polygons[[j]]@Polygons[[k]]@coords
								if (point.in.polygon(major_clade[i,"longitude"],major_clade[i,"latitude"],pol[,1],pol[,2]) == 1)
									{
										inState = TRUE
									}
							}
					}
				if (inState == TRUE) major_clade_boroughs[i,1] = "NY"
			}
	}
major_clade = cbind(major_clade, major_clade_boroughs)
major_clade = major_clade[which(!is.na(major_clade[,"location"])),c(1,3:5,2)]
fasta1 = scan(paste0(wd1,"/All_NY_sequences.fasta"), what="", sep="", quiet=T, blank.lines.skip=F)
sequenceIDs = gsub(">","",fasta1[which(grepl(">",fasta1))])
sequences = fasta1[which(!grepl(">",fasta1))]; fasta2 = c()
for (i in 1:dim(major_clade)[1])
	{
		major_clade_sequences[i] = sequences[which(sequenceIDs==major_clade[i,"trait"])]
		fasta2 = c(fasta2, paste0(">",major_clade[i,"trait"]), major_clade_sequences[i])
	}
write.table(major_clade, paste0(wd2,"/Major_clade_data.txt"), row.names=F, quote=F, sep="\t")
write(fasta2, paste0(wd2,"/Major_clade_seqs.fasta"))

tree = read.tree(paste0(wd1,"/",analysis,".tre"))
subtree = drop.tip(tree, tree$tip.label[!tree$tip.label%in%major_clade[,"trait"]])
mostRecentSamplingYear = max(decimal_date(ymd(major_clade[,"collection_date"])), na.rm=T)
node_height = max(nodeHeights(subtree)); tMRCA = mostRecentSamplingYear-node_height

		# Parameters of the BEAST analyses:
			# - substitution model: GTR+G
			# - TMRCA constrained to the age of the corresponding clade in the IQTREE+TreeTime tree
			# - fixed substitution rate (clock rate) estimated by Sam Hong (cfr. below)	
			# - coalescent model: skygrid (cut-off of 0.5 year, 50 grid points)
		
		# Analyses performed by Sam Hong to estimate the substitution rate to use:
			# - he randomly sampled maximum 3 unique sequences per day to get a total of 301 sequences (see "Screen_shot_Sam_1.png" for the TempEst root-to-tips regression)
			# - he then estimated the clock rate in BEAST with a GTR+G, strict clock and exponential growth models with the usual priors (--> clock rate of 8.413E-4)

	# B2. Preparing the integrated phylogeographic analyses on subsampled data sets

comparison = matrix(nrow=4, ncol=5); colnames(comparison) = gsub(" ","",NYboroughs@data[,"boro_name"])
row.names(comparison) = c("sampled_sequences","total_hospitalisations","ratio_samples_hospitalisations","sequences_to_subsample")
maxDate = as.character(max(dmy(gsub("\\/","-",metadata_NY_sequences[,"Clinical.Collection.Date"]))))
tab = read.csv("NY_epidemiological_data/NY_boroughs_COVID.csv", head=T)
hospitalisations = tab[which(tab[,"type"]=="ever-hospitalized"),]
hospitalisations = hospitalisations[which(grepl(maxDate,hospitalisations[,"timestamp"])),3:7]
colnames(hospitalisations) = gsub("_","",colnames(hospitalisations))
for (i in 1:dim(comparison)[2])
	{
		comparison[1,i] = sum(major_clade[,"location"]==colnames(comparison)[i], na.rm=T)
		comparison[2,i] = hospitalisations[1,which(colnames(hospitalisations)==tolower(colnames(comparison)[i]))]
		comparison[3,i] = comparison[1,i]/comparison[2,i]
	}
for (i in 1:dim(comparison)[2])
	{
		comparison[4,i] = round(comparison[3,which(comparison[3,]==min(comparison[3,]))]*comparison[2,i])
	}
for (n in 1:nberOfReplicates)
	{
		indices = c()
		for (i in 1:dim(comparison)[2])
			{
				buffer = which(major_clade[,"location"]==colnames(comparison)[i])
				indices = c(indices, sample(buffer, comparison["sequences_to_subsample",i], replace=F))
			}
		subsample1 = major_clade[as.vector(indices),]; subsample2 = sequences[as.vector(indices)]
		samplingDates = decimal_date(ymd(subsample1[,"collection_date"]))
		tree = read.tree(paste0(wd1,"/",analysis,".tre"))
		subtree = drop.tip(tree, tree$tip.label[!tree$tip.label%in%subsample1[,"trait"]])
		mostRecentSamplingYear = max(samplingDates, na.rm=T)
		node_height = max(nodeHeights(subtree)); tMRCA = mostRecentSamplingYear-node_height
		phylogeographicRuns = paste0(wd2,"/Replicate_DTA_analyses"); xml = c()
		write.tree(subtree, paste0(phylogeographicRuns,"/Replicate_",n,"_DTA.tre"))
		if (!file.exists(paste0(phylogeographicRuns,"/Replicate_",n,"_DTA.xml")))
			{
				template = scan(paste0(wd2,"/Template_file_DTA.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				template = gsub("TEMPLATE",paste0("Replicate_",n,"_DTA"), template)
				tree_txt = scan(paste0(phylogeographicRuns,"/Replicate_",n,"_DTA.tre"), what="", sep="\n", quiet=T)
				template = gsub("STARTING_TREE", tree_txt, template)
				template = gsub("NODE_HEIGHT", node_height, template)
				template = gsub("NODE_OFFSET", node_height-(7/365), template)
				sink(file=paste0(phylogeographicRuns,"/Replicate_",n,"_DTA.xml"))
				for (i in 1:length(template))
					{
						cat(template[i],"\n")
						if (grepl("Insert taxa block",template[i]))
							{
								cat(paste0("\t<taxa id=\"taxa\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat(paste0("\t\t<taxon id=\"",subsample1[j,"trait"],"\">","\n"))
										cat(paste0("\t\t\t<date value=\"",decimal_date(ymd(subsample1[j,"collection_date"])),"\" direction=\"forwards\" units=\"years\"/>","\n"))
										cat("\t\t\t<attr name=\"location\">\n")
										cat(paste0("\t\t\t\t",subsample1[j,"location"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t</taxon>\n")
									}
								cat("\t</taxa>","\n")
							}
						if (grepl("Insert clade block",template[i]))
							{
								cat(paste0("\t<taxa id=\"major_clade\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat(paste0("\t\t<taxon idref=\"",subsample1[j,"trait"],"\"/>","\n"))
									}
								cat("\t</taxa>","\n")
							}
						if (grepl("Insert alignment block",template[i]))
							{
								cat(paste0("\t<alignment id=\"alignment\" dataType=\"nucleotide\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat("\t\t<sequence>\n")
										cat(paste0("\t\t\t<taxon idref=\"",subsample1[j,"trait"],"\"/>","\n"))
										cat(paste0("\t\t\t",subsample2[j],"\n"))
										cat("\t\t</sequence>\n")
									}
								cat("\t</alignment>","\n")
							}
					}
				sink(NULL)
			}
		if (!file.exists(paste0(phylogeographicRuns,"/Replicate_",n,"_wTS.xml")))
			{
				template = scan(paste0(wd2,"/Template_file_wTS.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				template = gsub("TEMPLATE",paste0("Replicate_",n,"_wTS"), template)
				tree_txt = scan(paste0(phylogeographicRuns,"/Replicate_",n,"_DTA.tre"), what="", sep="\n", quiet=T)
				template = gsub("STARTING_TREE", tree_txt, template)
				template = gsub("NODE_HEIGHT", node_height, template)
				template = gsub("NODE_OFFSET", node_height-(7/365), template)
				sink(file=paste0(phylogeographicRuns,"/Replicate_",n,"_wTS.xml"))
				for (i in 1:length(template))
					{
						cat(template[i],"\n")
						if (grepl("Insert taxa block",template[i]))
							{
								cat(paste0("\t<taxa id=\"taxa\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat(paste0("\t\t<taxon id=\"",subsample1[j,"trait"],"\">","\n"))
										cat(paste0("\t\t\t<date value=\"",decimal_date(ymd(subsample1[j,"collection_date"])),"\" direction=\"forwards\" units=\"years\"/>","\n"))
										cat("\t\t\t<attr name=\"location\">\n")
										cat(paste0("\t\t\t\t",subsample1[j,"location"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t</taxon>\n")
									}
								cat("\t</taxa>","\n")
							}
						if (grepl("Insert clade block",template[i]))
							{
								cat(paste0("\t<taxa id=\"major_clade\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat(paste0("\t\t<taxon idref=\"",subsample1[j,"trait"],"\"/>","\n"))
									}
								cat("\t</taxa>","\n")
							}
						if (grepl("Insert alignment block",template[i]))
							{
								cat(paste0("\t<alignment id=\"alignment\" dataType=\"nucleotide\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat("\t\t<sequence>\n")
										cat(paste0("\t\t\t<taxon idref=\"",subsample1[j,"trait"],"\"/>","\n"))
										cat(paste0("\t\t\t",subsample2[j],"\n"))
										cat("\t\t</sequence>\n")
									}
								cat("\t</alignment>","\n")
							}
					}
				sink(NULL)
			}
		phylogeographicRuns = paste0(wd2,"/Replicate_RRW_analyses"); xml = c()
		write.tree(subtree, paste0(phylogeographicRuns,"/Replicate_",n,"_RRW.tre"))
		if (!file.exists(paste0(phylogeographicRuns,"/Replicate_",n,"_RRW.xml")))
			{
				template = scan(paste0(wd2,"/Template_file_RRW.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				template = gsub("TEMPLATE",paste0("Replicate_",n,"_RRW"), template)
				tree_txt = scan(paste0(phylogeographicRuns,"/Replicate_",n,"_RRW.tre"), what="", sep="\n", quiet=T)
				template = gsub("STARTING_TREE", tree_txt, template)
				template = gsub("NODE_HEIGHT", node_height, template)
				template = gsub("NODE_OFFSET", node_height-(7/365), template)
				sink(file=paste0(phylogeographicRuns,"/Replicate_",n,"_RRW.xml"))
				for (i in 1:length(template))
					{
						cat(template[i],"\n")
						if (grepl("Insert taxa block",template[i]))
							{
								cat(paste0("\t<taxa id=\"taxa\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat(paste0("\t\t<taxon id=\"",subsample1[j,"trait"],"\">","\n"))
										cat(paste0("\t\t\t<date value=\"",decimal_date(ymd(subsample1[j,"collection_date"])),"\" direction=\"forwards\" units=\"years\"/>","\n"))
										cat("\t\t\t<attr name=\"latitude\">\n")
										cat(paste0("\t\t\t\t",subsample1[j,"latitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"longitude\">\n")
										cat(paste0("\t\t\t\t",subsample1[j,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"location\">\n")
										cat(paste0("\t\t\t\t",subsample1[j,"latitude"]," ",subsample1[j,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t</taxon>\n")
									}
								cat("\t</taxa>","\n")
							}
						if (grepl("Insert clade block",template[i]))
							{
								cat(paste0("\t<taxa id=\"major_clade\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat(paste0("\t\t<taxon idref=\"",subsample1[j,"trait"],"\"/>","\n"))
									}
								cat("\t</taxa>","\n")
							}
						if (grepl("Insert alignment block",template[i]))
							{
								cat(paste0("\t<alignment id=\"alignment\" dataType=\"nucleotide\">","\n"))
								for (j in 1:dim(subsample1)[1])
									{
										cat("\t\t<sequence>\n")
										cat(paste0("\t\t\t<taxon idref=\"",subsample1[j,"trait"],"\"/>","\n"))
										cat(paste0("\t\t\t",subsample2[j],"\n"))
										cat("\t\t</sequence>\n")
									}
								cat("\t</alignment>","\n")
							}
					}
				sink(NULL)
			}
	}

	# B3. Analysing the output of the BSSVS and corresponding tip swapping analyses

wd = getwd()
setwd(paste0(wd,"/B_integrated_analyses/"))
if (!file.exists(paste0("MajorClade_1000.trees")))
	{
		trees1 = scan(paste0("MajorClade_DTA_1.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		indices1 = which(!grepl("tree STATE_", trees1)); indices2 = which(grepl("tree STATE_", trees1))
		interval = floor(length(indices2)/1000); indices2 = indices2[seq(interval,1000*interval,interval)]
		trees2 = c(trees1[c(indices1[1:(length(indices1)-1)],indices2)],"End;")
		write(trees2, paste0("MajorClade_1000.trees")) # then run TreeAnnotator
	}
if (!file.exists(paste0("MajorClade_DTA_4.log")))
	{
		log = scan(paste0("MajorClade_DTA_1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		interval = floor((length(log)-5)/1000); log = log[c(1:5,seq(6,6+(999*interval),interval))]
		write(log, paste0("MajorClade_DTA_4.log"))
	}
if (!file.exists(paste0("MajorClade_wTS_4.log")))
	{
		log = scan(paste0("MajorClade_wTS_1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		interval = floor((length(log)-5)/1000); log = log[c(1:5,seq(6,6+(999*interval),interval))]
		write(log, paste0("MajorClade_wTS_4.log"))
	}
if (!file.exists(paste0("MajorClade_DTA_5.log")))
	{
		log = scan(paste0("MajorClade_DTA_3.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		interval = floor((length(log)-4)/1000); log = log[c(1:4,seq(5,5+(999*interval),interval))]
		write(log, paste0("MajorClade_DTA_5.log"))
	}
if (!file.exists(paste0("MajorClade_wTS_5.log")))
	{
		log = scan(paste0("MajorClade_wTS_3.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		interval = floor((length(log)-4)/1000); log = log[c(1:4,seq(5,5+(999*interval),interval))]
		write(log, paste0("MajorClade_wTS_5.log"))
	}
if (!file.exists(paste0("MajorClade_DTA_5.txt")))
	{
		system(paste0("perl Getting_Markov_jumps.pl 0 < MajorClade_DTA_5.log > MajorClade_DTA_5.txt"))
		txt = scan(paste0("MajorClade_DTA_5.txt"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		write(txt[c(1,3:length(txt))], paste0("MajorClade_DTA_5.txt"))
	}
if (!file.exists(paste0("MajorClade_DTA_5.csv")))
	{
		tab = read.table(paste0("MajorClade_DTA_5.txt"), head=T)
		MJs = matrix(nrow=length(locations), ncol=length(locations))
		states = unique(tab[,"state"])
		for (i in 1:length(locations))
			{
				for (j in 1:length(locations))
					{
						sub = tab[which((tab[,"from"]==locations[i])&(tab[,"to"]==locations[j])),]
						MJs_post = c(); # print(c(i,j))
						for (k in 1:length(states))
							{
								MJ = which(sub[,"state"]==states[k])
								MJs_post = c(MJs_post, length(MJ))
							}
						MJs[i,j] = mean(MJs_post) # median ??
					}
			}
		row.names(MJs) = locations; colnames(MJs) = locations
		write.table(round(MJs,1), paste0("MajorClade_DTA_5.csv"), quote=F, sep=",")
	}	else		{
		MJs = read.csv(paste0("MajorClade_DTA_5.csv"), header=T)
	}
log1 = read.table(paste0("MajorClade_DTA_4.log"), header=T)
log2 = read.table(paste0("MajorClade_wTS_4.log"), header=T)
rates = matrix(nrow=length(locations), ncol=length(locations))
BFs1 = matrix(nrow=length(locations), ncol=length(locations))
BFs2 = matrix(nrow=length(locations), ncol=length(locations))
row.names(rates) = locations; colnames(rates) = locations
row.names(BFs1) = locations; colnames(BFs1) = locations
row.names(BFs2) = locations; colnames(BFs2) = locations
for (i in 1:length(locations))
	{
		for (j in 1:length(locations))
			{
				if (i != j)
					{
						colName = paste0("location.indicators.",locations[i],".",locations[j])
						index1 = which(colnames(log1)==colName); index2 = which(colnames(log2)==colName)
						p = sum(log1[,index1]==1)/dim(log1)[1]
						K = 20 # length(locations)*(length(locations)-1) # K shoulf be divided by 2 if "symetric" case
						q = (log(2)+K-1)/(K*(K-1))
						BFs1[i,j] = (p/(1-p))/(q/(1-q))
						p1 = sum(log1[,index1]==1)/dim(log1)[1]
						p2 = sum(log2[,index2]==1)/dim(log2)[1]
						BFs2[i,j] = (p1/(1-p1))/(p2/(1-p2))
						index3 = which(colnames(log1)==paste0("location.rates.",locations[i],".",locations[j]))
						rates[i,j] = median(log1[,index3]*log1[,index1])
					}
			}
	}
MJs1 = MJs; MJs2 = MJs; rates[BFs2<3] = NA; MJs1[BFs1<3] = NA; MJs2[BFs2<3] = NA; BFs1 = round(BFs1, 1); BFs2 = round(BFs2, 1)
write.table(BFs1, paste0("MajorClade_BFs_1.csv"), sep=",", quote=F)
write.table(BFs2, paste0("MajorClade_BFs_2.csv"), sep=",", quote=F)
write.table(rates, paste0("MajorClade_rates.csv"), sep=",", quote=F)
write.table(MJs2, paste0("MajorClade_MJs_1.csv"), sep=",", quote=F)
write.table(MJs2, paste0("MajorClade_MJs_2.csv"), sep=",", quote=F)

wd = getwd()
setwd(paste0(wd,"/B_integrated_analyses/Replicate_DTA_analyses/"))
locations = c("Bronx","Brooklyn","Manhattan","Queens","StatenIsland")
burnIns_DTA = c(1001,301,1041,301,301,241,301,1001,201,4401)
burnIns_wTS = c(2001,1501,1801,1301,261,1701,201,1301,801,2201)
for (n in 1:nberOfReplicates)
	{
		if (!file.exists(paste0("Replicate_",n,"_1000.trees")))
			{
				trees1 = scan(paste0("Replicate_",n,"_DTA_1.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				indices1 = which(!grepl("tree STATE_", trees1)); indices2 = which(grepl("tree STATE_", trees1))
				indices2 = indices2[(burnIns_DTA[n]+1):length(indices2)]; interval = floor(length(indices2)/1000)
				indices2 = indices2[seq(interval,1000*interval,interval)]
				trees2 = c(trees1[c(indices1[1:(length(indices1)-1)],indices2)],"End;")
				write(trees2, paste0("Replicate_",n,"_1000.trees"))
			}
	}
for (n in 1:nberOfReplicates)
	{
		system(paste0("BEAST_1104_program/bin/treeannotator -burninTrees 0 -heights keep Replicate_",n,"_1000.trees Replicate_",n,"_MCC.tree"), ignore.stdout=F, ignore.stderr=F)
	}
for (n in 1:nberOfReplicates)
	{
		if (!file.exists(paste0("Replicate_",n,"_DTA_4.log")))
			{
				log = scan(paste0("Replicate_",n,"_DTA_1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				log = log[c(1:5,(burnIns_DTA[n]+1):length(log))]; interval = floor((length(log)-5)/1000)
				log = log[c(1:5,seq(6,6+(999*interval),interval))]; write(log, paste0("Replicate_",n,"_DTA_4.log"))
			}
		if (!file.exists(paste0("Replicate_",n,"_wTS_4.log")))
			{
				log = scan(paste0("Replicate_",n,"_wTS_1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				log = log[c(1:5,(burnIns_DTA[n]+1):length(log))]; interval = floor((length(log)-5)/1000)
				log = log[c(1:5,seq(6,6+(999*interval),interval))]; write(log, paste0("Replicate_",n,"_wTS_4.log"))
			}
		if (!file.exists(paste0("Replicate_",n,"_DTA_5.log")))
			{
				log = scan(paste0("Replicate_",n,"_DTA_3.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				log = log[c(1:4,(burnIns_DTA[n]+1):length(log))]; interval = floor((length(log)-4)/1000)
				log = log[c(1:4,seq(5,5+(999*interval),interval))]; write(log, paste0("Replicate_",n,"_DTA_5.log"))
			}
		if (!file.exists(paste0("Replicate_",n,"_wTS_5.log")))
			{
				log = scan(paste0("Replicate_",n,"_wTS_3.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				log = log[c(1:4,(burnIns_DTA[n]+1):length(log))]; interval = floor((length(log)-4)/1000)
				log = log[c(1:4,seq(5,5+(999*interval),interval))]; write(log, paste0("Replicate_",n,"_wTS_5.log"))
			}
	}
for (n in 1:nberOfReplicates)
	{
		if (!file.exists(paste0("Replicate_",n,"_DTA_5.txt")))
			{
				system(paste0("perl Getting_Markov_jumps.pl 0 < Replicate_",n,"_DTA_5.log > Replicate_",n,"_DTA_5.txt"))
				txt = scan(paste0("Replicate_",n,"_DTA_5.txt"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				write(txt[c(1,3:length(txt))], paste0("Replicate_",n,"_DTA_5.txt"))
			}
		if (!file.exists(paste0("Replicate_",n,"_DTA_5.csv")))
			{
				tab = read.table(paste0("Replicate_",n,"_DTA_5.txt"), head=T)
				MJs = matrix(nrow=length(locations), ncol=length(locations))
				states = unique(tab[,"state"])
				for (i in 1:length(locations))
					{
						for (j in 1:length(locations))
							{
								sub = tab[which((tab[,"from"]==locations[i])&(tab[,"to"]==locations[j])),]
								MJs_post = c(); # print(c(i,j))
								for (k in 1:length(states))
									{
										MJ = which(sub[,"state"]==states[k])
										MJs_post = c(MJs_post, length(MJ))
									}
								MJs[i,j] = mean(MJs_post) # median ??
							}
					}
				row.names(MJs) = locations; colnames(MJs) = locations
				write.table(round(MJs,1), paste0("Replicate_",n,"_DTA_5.csv"), quote=F, sep=",")
			}	else		{
				MJs = read.csv(paste0("Replicate_",n,"_DTA_5.csv"), header=T)
			}
		log1 = read.table(paste0("Replicate_",n,"_DTA_4.log"), header=T)
		log2 = read.table(paste0("Replicate_",n,"_wTS_4.log"), header=T)
		rates = matrix(nrow=length(locations), ncol=length(locations))
		BFs1 = matrix(nrow=length(locations), ncol=length(locations))
		BFs2 = matrix(nrow=length(locations), ncol=length(locations))
		row.names(rates) = locations; colnames(rates) = locations
		row.names(BFs1) = locations; colnames(BFs1) = locations
		row.names(BFs2) = locations; colnames(BFs2) = locations
		for (i in 1:length(locations))
			{
				for (j in 1:length(locations))
					{
						if (i != j)
							{
								colName = paste0("location.indicators.",locations[i],".",locations[j])
								index1 = which(colnames(log1)==colName); index2 = which(colnames(log2)==colName)
								p = sum(log1[,index1]==1)/dim(log1)[1]
								K = 20 # length(locations)*(length(locations)-1) # K shoulf be divided by 2 if "symetric" case
								q = (log(2)+K-1)/(K*(K-1))
								BFs1[i,j] = (p/(1-p))/(q/(1-q))
								p1 = sum(log1[,index1]==1)/dim(log1)[1]
								p2 = sum(log2[,index2]==1)/dim(log2)[1]
								BFs2[i,j] = (p1/(1-p1))/(p2/(1-p2))
								index3 = which(colnames(log1)==paste0("location.rates.",locations[i],".",locations[j]))
								rates[i,j] = median(log1[,index3]*log1[,index1])
							}
					}
			}
		MJs1 = MJs; MJs2 = MJs; rates[BFs2<3] = NA; MJs1[BFs1<3] = NA; MJs2[BFs2<3] = NA; BFs1 = round(BFs1, 1); BFs2 = round(BFs2, 1)
		write.table(BFs1, paste0("Replicate_",n,"_BFs_1.csv"), sep=",", quote=F)
		write.table(BFs2, paste0("Replicate_",n,"_BFs_2.csv"), sep=",", quote=F)
		write.table(rates, paste0("Replicate_",n,"_rates.csv"), sep=",", quote=F)
		write.table(MJs1, paste0("Replicate_",n,"_MJs_1.csv"), sep=",", quote=F)
		write.table(MJs2, paste0("Replicate_",n,"_MJs_2.csv"), sep=",", quote=F)
	}
setwd(wd)

if (showingPlots)
	{
		MJs = read.csv(paste0("B_integrated_analyses/Replicate_DTA_analyses/Replicate_",1,"_MJs_1.csv"), head=T)
		mat = as.matrix(MJs); diag(mat) = NA; minVals = min(mat, na.rm=T); maxVals = max(mat, na.rm=T)
		for (n in 2:nberOfReplicates)
			{
				MJs = read.csv(paste0("B_integrated_analyses/Replicate_DTA_analyses/Replicate_",n,"_MJs_1.csv"), head=T)
				mat = as.matrix(MJs); diag(mat) = NA; minVals = min(mat, na.rm=T); maxVals = max(mat, na.rm=T)
				if (minVals > min(mat,na.rm=T)) minVals = min(mat, na.rm=T)
				if (maxVals < max(mat,na.rm=T)) maxVals = max(mat, na.rm=T)
			}
		BayesFactors = 1; BayesFactors = 2
		if (BayesFactors == 1) pdf(paste0("Figure_E1_NEW.pdf"), width=11.0, height=4.5) # dev.new(width=11.0, height=4.5)
		if (BayesFactors == 2) pdf(paste0("Figure_E2_NEW.pdf"), width=11.0, height=4.5) # dev.new(width=11.0, height=4.5)
		par(mfrow=c(2,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		for (n in 1:nberOfReplicates)
			{
				mat = as.matrix(read.csv(paste0("B_integrated_analyses/Replicate_DTA_analyses/Replicate_",n,"_MJs_",BayesFactors,".csv"), head=T))
				colNames = gsub(" ","",NYboroughs@data[,"boro_name"]); mat = mat[colNames,colNames]
				plot(NYboroughs, col=col_boroughs, border="gray60", lwd=0.5)
				points(centroids, cex=10*((diag(mat)-minVals)/(maxVals1-minVals)), pch=16, col="#4D4D4D50")
				points(centroids, cex=10*((vec-minVals)/(maxVals1-minVals)), pch=1, col="#4D4D4D75", lwd=0.5, lty=2)
				for (i in 1:dim(NYboroughs)[1])
					{
						for (j in 1:dim(NYboroughs)[1])
							{
								if ((!is.na(mat[i,j]))&&(i != j)&&(mat[i,j]>=1))
									{
										LWD = (((mat[i,j]-minVals)/(maxVals-minVals))*2)+0.1; arrow = (0.1*(mat[i,j]/maxVals))+0.04
										curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
													lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
									}
							}
					}
				if (n == 1)
					{
						vS = 5; LWD = (((vS-minVals)/(maxVals-minVals))*2)+0.1; arrow = (0.1*(vS/maxVals))+0.04
						curvedarrow(cbind(-74.20,40.900), cbind(-74.12,40.900), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 10; LWD = (((vS-minVals)/(maxVals-minVals))*2)+0.1; arrow = (0.1*(vS/maxVals))+0.04
						curvedarrow(cbind(-74.20,40.875), cbind(-74.12,40.875), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 20; LWD = (((vS-minVals)/(maxVals-minVals))*2)+0.1; arrow = (0.1*(vS/maxVals))+0.04
						curvedarrow(cbind(-74.20,40.850), cbind(-74.12,40.850), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						mtext("5", at=-74.11, line=-1.90, adj=0, col="gray30", cex=0.7)
						mtext("10", at=-74.11, line=-2.80, adj=0, col="gray30", cex=0.7)
						mtext("20", at=-74.11, line=-3.70, adj=0, col="gray30", cex=0.7)
					}
			}
		dev.off()
	}

	# B4. Extracting spatio-temporal information embedded in posterior trees (DTA approach)

wd = getwd()
setwd(paste0(wd,"/B_integrated_analyses/Replicate_DTA_analyses/"))
for (n in 1:nberOfReplicates)
	{
		dir.create(file.path(paste0("Replicate_",n,"_DTA_ext")), showWarnings=F)
		if (!file.exists(paste0("Replicate_",n,"_DTA_ext/TreeExtractions_1.csv")))
			{
				trees = readAnnotatedNexus(paste0("Replicate_",n,"_1000.trees"))
				for (i in 1:length(trees))
					{
						tree = trees[[i]]
						tab = matrix(nrow=dim(tree$edge)[1], ncol=4)
						colnames(tab) = c("node1","node2","startLoc","endLoc")
						tab[,"node1"] = tree$edge[,1]; tab[,"node2"] = tree$edge[,2]
						for (j in 1:dim(tree$edge)[1])
							{
								tab[j,"endLoc"] = tree$annotations[[j]]$location
								index = which(tree$edge[,2]==tree$edge[j,1])
								if (length(index) == 1)
									{
										tab[j,"startLoc"] = tree$annotations[[index]]$location
									}	else		{
										if (!tree$edge[j,1]%in%tree$edge[,2])
											{
												tab[j,"startLoc"] = tree$root.annotation$location
											}
									}
							}
						write.csv(tab, paste0("Replicate_",n,"_DTA_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
					}
			}
	}
setwd(wd)

	# B5. Extracting spatio-temporal information embedded in MCC/posterior trees (RRW approach)

source("MCC_tree_extraction.r")
metadata = read.table("B_integrated_analyses/Major_clade_data.txt", head=T)
mostRecentSamplingDates = rep(NA, nberOfReplicates); wd = getwd()
setwd(paste0(wd,"/B_integrated_analyses/Replicate_RRW_analyses/"))
burnIns = rep(NA, nberOfReplicates)
for (n in 1:nberOfReplicates)
	{
		allTrees = scan(file=paste0("Replicate_",n,"_RRW.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		burnIns[n] = round(sum(grepl("tree STATE_",allTrees))/10)+1
	}
for (n in 1:nberOfReplicates)
	{
		if (!file.exists(paste0("Replicate_",n,"_MCC.tree")))
			{
				system(paste0("BEAST_1104_program/bin/treeannotator -burninTrees ",burnIns[n]," -heights keep Replicate_",n,"_RRW.trees Replicate_",n,"_MCC.tree"), ignore.stdout=F, ignore.stderr=F)
			}
	}
for (n in 1:nberOfReplicates)
	{
		if (!file.exists(paste0("Replicate_",n,"_MCC.csv")))
			{
				mcc_tre = readAnnotatedNexus(paste0("Replicate_",n,"_MCC.tree"))
				indices = which(metadata[,"trait"]%in%mcc_tre$tip.label)
				if (length(indices) != length(mcc_tre$tip.label))
					{
						print(n)
					}	else	{
						mostRecentSamplingDatum = max(decimal_date(ymd(metadata[indices,"collection_date"])))
						mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum); mostRecentSamplingDates[n] = mostRecentSamplingDatum
						write.csv(mcc_tab, paste0("Replicate_",n,"_MCC.csv"), row.names=F, quote=F)
					}
			}
	}
randomSampling = FALSE; nberOfTreesToSample = 1000; coordinateAttributeName = "location"; nberOfCores = 10
for (n in 1:nberOfReplicates)
	{
		localTreesDirectory = paste0("Replicate_",n,"_RRW_ext")
		if (!file.exists(paste0(localTreesDirectory,"/TreeExtractions_1.csv")))
			{
				burnIn = burnIns[n]; mostRecentSamplingDatum = mostRecentSamplingDates[n]
				allTrees = scan(file=paste0("Replicate_",n,"_RRW.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
			}
	}
for (n in 1:nberOfReplicates)
	{
		mcc_tab = read.csv(paste0("Replicate_",n,"_MCC.csv"), head=T)
		if (!"tipLabel"%in%colnames(mcc_tab))
			{
				tipLabels = matrix(nrow=dim(mcc_tab)[1], ncol=1); colnames(tipLabels) = "tipLabel"
				for (j in 1:dim(mcc_tab)[1])
					{
						if (!mcc_tab[j,"node2"]%in%mcc_tab[,"node1"])
							{
								index = which((round(metadata[,"longitude"],5)==round(mcc_tab[j,"endLon"],5))&(round(metadata[,"latitude"],5)==round(mcc_tab[j,"endLat"],5)))
								if (length(index) != 1)
									{
										print(c(i,j))
									}	else	{
										tipLabels[j,1] = metadata[index,"trait"]
									}
							}	
					}
				mcc_tab = cbind(mcc_tab, tipLabels)
				write.csv(mcc_tab, paste0("Replicate_",n,"_MCC.csv"), row.names=F, quote=F)
			}
	}
for (n in 1:nberOfReplicates)
	{
		for (j in 1:nberOfTreesToSample)
			{
				tab = read.csv(paste0("Replicate_",n,"_RRW_ext/TreeExtractions_",j,".csv"), head=T)
				if (!"tipLabel"%in%colnames(tab))
					{
						tipLabels = matrix(nrow=dim(tab)[1], ncol=1); colnames(tipLabels) = "tipLabel"
						for (k in 1:dim(tab)[1])
							{
								if (!tab[k,"node2"]%in%tab[,"node1"])
									{
										index = which((round(metadata[,"longitude"],5)==round(tab[k,"endLon"],5))&(round(metadata[,"latitude"],5)==round(tab[k,"endLat"],5)))
										if (length(index) != 1)
											{
												print(c(i,j,k))
											}	else	{
												tipLabels[k,1] = metadata[index,"trait"]
											}
									}	
							}
						tab = cbind(tab, tipLabels)
						write.csv(tab, paste0("Replicate_",n,"_RRW_ext/TreeExtractions_",j,".csv"), row.names=F, quote=F)
					}
			}
	}
setwd(wd)

	# B6. Generating dispersal history graphs (RRW approach, MCC trees and 80% HPD polygons)

if (showingPlots)
	{
		wd = getwd(); nberOfReplicates = 10; mccs = list(); polygons_list = list()
		setwd(paste0(wd,"/B_integrated_analyses/Replicate_RRW_analyses/"))
		for (n in 1:nberOfReplicates)
			{
				mccs[[n]] = read.csv(paste0("Replicate_",n,"_MCC.csv"), head=T)
				localTreesDirectory = paste0("Replicate_",n,"_RRW_ext"); startDatum = min(mccs[[n]][,"startYear"])
				nberOfExtractionFiles = 1000; percentage = 80; prob = percentage/100; precision = 1/(365/7)
				polygons_list[[n]] = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
			}
		minYear = min(mccs[[1]][,"startYear"]); maxYear = max(mccs[[1]][,"endYear"])
		for (n in 2:nberOfReplicates)
			{
				if (minYear > min(mccs[[n]][,"startYear"])) minYear = min(mccs[[n]][,"startYear"])
				if (maxYear < max(mccs[[n]][,"endYear"])) maxYear = max(mccs[[n]][,"endYear"])
			}
		colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(141)[16:116])
		setwd(wd); onlyInternalNodesOfTipBranches = FALSE
		if (onlyInternalNodesOfTipBranches != TRUE) pdf(paste0("Figure_F_NEW1.pdf"), width=11.0, height=4.5)
		if (onlyInternalNodesOfTipBranches == TRUE) pdf(paste0("Figure_F_NEW2.pdf"), width=11.0, height=4.5)
		par(mfrow=c(2,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); cexNode = 0.7; LWD = 1.0; croppingPolygons = TRUE
		for (n in 1:nberOfReplicates)
			{
				polygons = polygons_list[[n]]; mcc = mccs[[n]]; selectedBranches = 1:dim(mcc)[1]
				startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
				startYears_colours = colourScale[startYears_indices]
				endYears_colours = colourScale[endYears_indices]
				polygons_colours = rep(NA, length(polygons_list[[]]))
				for (i in 1:length(polygons))
					{
						date = as.numeric(names(polygons[[i]])); polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
						polygons_colours[i] = paste0(colourScale[polygon_index],"40")
					}
				plot(NYboroughs, col="gray90", border="white", lwd=0.5)
				for (i in length(polygons):1)
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]
						if (croppingPolygons == TRUE)
							{
								pol = crop(pol, NYboroughs)
							}
						plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
					}
				plot(NYboroughs, col=NA, border="white", lwd=0.5, add=T)
				for (i in selectedBranches)
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
								    arr.width=0, lwd=0.1, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
				for (i in rev(selectedBranches))
					{
						if (onlyInternalNodesOfTipBranches != TRUE)
							{
								if (!mcc[i,"node1"]%in%mcc[selectedBranches,"node2"])
									{
										points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
										points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=cexNode, lwd=0.2)
									}
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=cexNode, lwd=0.2)
							}	else	{
								if (!mcc[i,"node2"]%in%mcc[selectedBranches,"node1"])
									{
										points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
										points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=cexNode, lwd=0.2)
									}
							}			
					}
				if (n == nberOfReplicates)
					{
						selectedDates = decimal_date(ymd(c("2020-03-01","2020-04-01","2020-05-01")))
						selectedLabels = c("01-03","01-04","01-05")	
						rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.970,0.990,0.1,0.9),
							 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=F,
				  			 axis.args=list(cex.axis=0.8, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.5,0),
				   			 at=selectedDates, labels=selectedLabels))
				     }
			}
		dev.off()
	}

	# B7. Analysing the lineage dispersal events among NY boroughs (DTA approach)

nberOfExtractionFiles = 1000
if (!file.exists("B_integrated_analyses/Replicate_DTA_analyses/Pairwise_matrices.rds"))
	{
		matrices_list = list()
		for (n in 1:nberOfReplicates)
			{
				matrices = list()
				for (i in 1:nberOfExtractionFiles)
					{
						mat = matrix(0, nrow=dim(NYboroughs@data)[1], ncol=dim(NYboroughs@data)[1])
						tab = read.csv(paste0("B_integrated_analyses/Replicate_DTA_analyses/Replicate_",n,"_DTA_ext/TreeExtractions_",i,".csv"), head=T)
						for (j in 1:dim(tab)[1])
							{
								index1 = which(gsub(" ","",NYboroughs@data[,"boro_name"])==tab[j,"startLoc"])
								index2 = which(gsub(" ","",NYboroughs@data[,"boro_name"])==tab[j,"endLoc"])
								mat[index1,index2] = mat[index1,index2]+1
							}
						matrices[[i]] = mat
					}
				matrices_list[[n]] = matrices
			}
		saveRDS(matrices_list, "B_integrated_analyses/Replicate_DTA_analyses/Pairwise_matrices.rds")
	}	else		{
		matrices_list = readRDS("B_integrated_analyses/Replicate_DTA_analyses/Pairwise_matrices.rds")
	}
matrices_mean = list()
for (n in 1:nberOfReplicates)
	{
		mat = matrix(0, nrow=dim(NYboroughs@data)[1], ncol=dim(NYboroughs@data)[1])
		for (i in 1:nberOfExtractionFiles) mat = mat+matrices_list[[n]][[i]]
		matrices_mean[[n]] = mat/nberOfExtractionFiles
	}
if (showingPlots)
	{
		centroids = coordinates(NYboroughs)	
		pdf(paste0("Figure_G_NEW.pdf"), width=11.0, height=4.5) # dev.new(width=11.0, height=3.9)
		par(mfrow=c(2,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		minVals1 = min(diag(matrices_mean[[1]])); maxVals1 = max(diag(matrices_mean[[1]]))
		mat = matrices_mean[[1]]; diag(mat) = NA; minVals2 = min(mat, na.rm=T); maxVals2 = max(mat, na.rm=T)
		for (n in 2:nberOfReplicates)
			{
				mat1 = matrices_mean[[n]]; mat2 = mat1; diag(mat2) = NA
				if (minVals1 > min(diag(mat1))) minVals1 = min(diag(mat1))
				if (maxVals1 < max(diag(mat1))) maxVals1 = max(diag(mat1))
				if (minVals2 > min(mat2,na.rm=T)) minVals2 = min(mat2,na.rm=T)
				if (maxVals2 < max(mat2,na.rm=T)) maxVals2 = max(mat2,na.rm=T)
			}
		for (n in 1:nberOfReplicates)
			{
				mat = matrices_mean[[n]]
				plot(NYboroughs, col=col_boroughs, border="gray60", lwd=0.5)
				points(centroids, cex=10*((diag(mat)-minVals1)/(maxVals1-minVals1)), pch=16, col="#4D4D4D50")
				for (i in 1:dim(NYboroughs)[1])
					{
						for (j in 1:dim(NYboroughs)[1])
							{
								if ((i != j)&(mat[i,j]>=1))
									{
										LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(mat[i,j]/maxVals2))+0.04
										curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
													lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
									}
							}
					}
				if (n == 1)
					{
						vS = 5; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.900), cbind(-74.12,40.900), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 10; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.875), cbind(-74.12,40.875), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 20; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.850), cbind(-74.12,40.850), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						mtext("5", at=-74.11, line=-1.90, adj=0, col="gray30", cex=0.7)
						mtext("10", at=-74.11, line=-2.80, adj=0, col="gray30", cex=0.7)
						mtext("20", at=-74.11, line=-3.70, adj=0, col="gray30", cex=0.7)
						points(cbind(rep(-74.15,4),rep(40.75,4)), cex=10*((seq(40,120,40)-minVals1)/(maxVals1-minVals1)), pch=1, col="#4D4D4D", lwd=0.3)
					}
			}		
		dev.off()
	}

	# B8. Analysing the lineage dispersal events among NY boroughs (RRW approach)

nberOfExtractionFiles = 1000
if (!file.exists("B_integrated_analyses/Replicate_RRW_analyses/Pairwise_matrices_1.rds"))
	{
		matrices_list1 = list() # considering all the branches
		matrices_list2 = list() # only considering tip branches
		for (n in 1:nberOfReplicates)
			{
				matrices1 = list(); matrices2 = list()
				for (i in 1:nberOfExtractionFiles)
					{
						mat1 = matrix(0, nrow=dim(NYboroughs@data)[1], ncol=dim(NYboroughs@data)[1])
						mat2 = matrix(0, nrow=dim(NYboroughs@data)[1], ncol=dim(NYboroughs@data)[1])
						tab = read.csv(paste0("B_integrated_analyses/Replicate_RRW_analyses/Replicate_",n,"_RRW_ext/TreeExtractions_",i,".csv"), head=T)
						for (j in 1:dim(tab)[1])
							{
								polIndex1 = c(); polIndex2 = c()
								for (k in 1:length(NYboroughs@polygons))
									{
										for (l in 1:length(NYboroughs@polygons[[k]]@Polygons))
											{
												pol = NYboroughs@polygons[[k]]@Polygons[[l]]
												if (point.in.polygon(tab[j,"startLon"],tab[j,"startLat"],pol@coords[,1],pol@coords[,2]) == 1)
													{
														polIndex1 = c(polIndex1, k)
													}
												if (point.in.polygon(tab[j,"endLon"],tab[j,"endLat"],pol@coords[,1],pol@coords[,2]) == 1)
													{
														polIndex2 = c(polIndex2, k)
													}
											}
									}
								if ((length(polIndex1)==1)&(length(polIndex2)==1))
									{
										mat1[polIndex1,polIndex2] = mat1[polIndex1,polIndex2] + 1
										if (!tab[j,"node2"]%in%tab[,"node1"])
											{
												mat2[polIndex1,polIndex2] = mat2[polIndex1,polIndex2] + 1
											}
									}
							}
						matrices1[[i]] = mat1; matrices2[[i]] = mat2
					}
				matrices_list1[[n]] = matrices1; matrices_list2[[n]] = matrices2
			}
		saveRDS(matrices_list1, "B_integrated_analyses/Replicate_RRW_analyses/Pairwise_matrices_1.rds")
		saveRDS(matrices_list2, "B_integrated_analyses/Replicate_RRW_analyses/Pairwise_matrices_2.rds")
	}	else		{
		matrices_list1 = readRDS("B_integrated_analyses/Replicate_RRW_analyses/Pairwise_matrices_1.rds")
		matrices_list2 = readRDS("B_integrated_analyses/Replicate_RRW_analyses/Pairwise_matrices_2.rds")
	}
if (!file.exists("B_integrated_analyses/Replicate_RRW_analyses/Borough_branches_1.rds"))
	{
		branches_list1 = list() # number of branches within clades connecting sequences sampled within a unique borough
		branches_list2 = list() # number of tip branches connecting sequences sampled within a unique borough
		for (n in 1:nberOfReplicates)
			{
				branches1 = list(); branches2 = list()
				for (i in 1:nberOfExtractionFiles)
					{
						vals1 = rep(NA, dim(NYboroughs@data)[1]); vals2 = rep(NA, dim(NYboroughs@data)[1])
						tab = read.csv(paste0("B_integrated_analyses/Replicate_RRW_analyses/Replicate_",n,"_RRW_ext/TreeExtractions_",i,".csv"), head=T)
						for (j in 1:length(NYboroughs@polygons))
							{
								branches = c()
								for (k in 1:dim(tab)[1])
									{
										if (!tab[k,"node2"]%in%tab[,"node1"])
											{
												polIndex = c()
												for (l in 1:length(NYboroughs@polygons[[j]]@Polygons))
													{
														pol = NYboroughs@polygons[[j]]@Polygons[[l]]
														if (point.in.polygon(tab[k,"endLon"],tab[k,"endLat"],pol@coords[,1],pol@coords[,2]) == 1) polIndex = l
													}
												if (length(polIndex) == 1) branches = c(branches,k)
											}
									}
								sub = tab[branches,]	
								nodesToExplore = tab[which(tab[,"node2"]%in%sub[,"node1"]),"node2"]
								buffer1 = nodesToExplore			
								while (length(buffer1) != 0)
									{
										buffer1 = c(); buffer2 = c()
										for (k in 1:length(nodesToExplore))
											{
												if (sum(sub[,"node1"]==nodesToExplore[k]) == 2)
													{
														sub = rbind(sub, tab[which(tab["node2"]==nodesToExplore[k]),])
														buffer1 = c(buffer1, tab[which(tab["node2"]==nodesToExplore[k]),"node1"])
													}	else		{
														buffer2 = c(buffer2, nodesToExplore[k])
													}	
											}
										nodesToExplore = c(buffer1, buffer2)
									}
								vals1[j] = dim(sub)[1]
								sub = sub[which(!sub[,"node2"]%in%sub[,"node1"]),]; buffer = c()
								for (k in 1:dim(sub)[1])
									{
										if (sum(sub[,"node1"]==sub[k,"node1"]) == 2) buffer = c(buffer, k)
									}
								vals2[j] = length(buffer)
							}
						branches1[[i]] = vals1; branches2[[i]] = vals2
					}
				branches_list1[[n]] = branches1; branches_list2[[n]] = branches2
			}
		saveRDS(branches_list1, "B_integrated_analyses/Replicate_RRW_analyses/Borough_branches_1.rds")
		saveRDS(branches_list2, "B_integrated_analyses/Replicate_RRW_analyses/Borough_branches_2.rds")
	}	else		{
		branches_list1 = readRDS("B_integrated_analyses/Replicate_RRW_analyses/Borough_branches_1.rds")
		branches_list2 = readRDS("B_integrated_analyses/Replicate_RRW_analyses/Borough_branches_2.rds")
	}
matrices_mean1 = list(); branches_mean1 = list() # considering all the branches
matrices_mean2 = list(); branches_mean2 = list() # only considering tip branches
for (n in 1:nberOfReplicates)
	{
		mat1 = matrix(0, nrow=dim(NYboroughs@data)[1], ncol=dim(NYboroughs@data)[1])
		mat2 = matrix(0, nrow=dim(NYboroughs@data)[1], ncol=dim(NYboroughs@data)[1])
		for (i in 1:nberOfExtractionFiles)
			{
				mat1 = mat1+matrices_list1[[n]][[i]]; mat2 = mat2+matrices_list2[[n]][[i]]
			}
		matrices_mean1[[n]] = mat1/nberOfExtractionFiles
		matrices_mean2[[n]] = mat2/nberOfExtractionFiles
		vec1 = rep(0, dim(NYboroughs@data)[1])
		vec2 = rep(0, dim(NYboroughs@data)[1])
		for (i in 1:nberOfExtractionFiles)
			{
				vec1 = vec1+branches_list1[[n]][[i]]; vec2 = vec2+branches_list2[[n]][[i]]
			}
		branches_mean1[[n]] = vec1/nberOfExtractionFiles
		branches_mean2[[n]] = vec2/nberOfExtractionFiles
	}
if (showingPlots)
	{
		centroids = coordinates(NYboroughs); plottingBoroughCladeBranches = FALSE
		pdf(paste0("Figure_H_NEW1.pdf"), width=11.0, height=4.5) # dev.new(width=11.0, height=3.9)
		par(mfrow=c(2,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		minVals1 = min(diag(matrices_mean1[[1]])); maxVals1 = max(diag(matrices_mean1[[1]]))
		mat = matrices_mean1[[1]]; diag(mat) = NA; minVals2 = min(mat, na.rm=T); maxVals2 = max(mat, na.rm=T)
		vec = branches_mean1[[1]]; minVals3 = min(vec, na.rm=T); maxVals3 = max(vec, na.rm=T)
		for (n in 2:nberOfReplicates)
			{
				mat1 = matrices_mean1[[n]]; mat2 = mat1
				diag(mat2) = NA; mat3 = branches_mean1[[n]]
				if (minVals1 > min(diag(mat1))) minVals1 = min(diag(mat1))
				if (maxVals1 < max(diag(mat1))) maxVals1 = max(diag(mat1))
				if (minVals2 > min(mat2,na.rm=T)) minVals2 = min(mat2,na.rm=T)
				if (maxVals2 < max(mat2,na.rm=T)) maxVals2 = max(mat2,na.rm=T)
				if (minVals3 > min(mat3,na.rm=T)) minVals3 = min(mat3,na.rm=T)
				if (maxVals3 < max(mat3,na.rm=T)) maxVals3 = max(mat3,na.rm=T)
			}
		for (n in 1:nberOfReplicates)
			{
				mat = matrices_mean1[[n]]; vec = branches_mean1[[n]]
				plot(NYboroughs, col=col_boroughs, border="gray60", lwd=0.5)
				points(centroids, cex=10*((diag(mat)-minVals1)/(maxVals1-minVals1)), pch=16, col="#4D4D4D50")
				if (plottingBoroughCladeBranches == TRUE)
					{
						points(centroids, cex=10*((vec-minVals1)/(maxVals1-minVals1)), pch=1, col="#4D4D4D75", lwd=0.5, lty=2)
					}
				for (i in 1:dim(NYboroughs)[1])
					{
						for (j in 1:dim(NYboroughs)[1])
							{
								if ((i != j)&(mat[i,j]>=1))
									{
										LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(mat[i,j]/maxVals2))+0.04
										curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
													lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
									}
							}
					}
				if (n == 1)
					{
						vS = 5; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.900), cbind(-74.12,40.900), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 10; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.875), cbind(-74.12,40.875), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 20; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.850), cbind(-74.12,40.850), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						mtext("5", at=-74.11, line=-1.90, adj=0, col="gray30", cex=0.7)
						mtext("10", at=-74.11, line=-2.80, adj=0, col="gray30", cex=0.7)
						mtext("20", at=-74.11, line=-3.70, adj=0, col="gray30", cex=0.7)
						points(cbind(rep(-74.15,4),rep(40.75,4)), cex=10*((seq(20,80,20)-minVals1)/(maxVals1-minVals1)), pch=1, col="#4D4D4D", lwd=0.3)
					}
			}		
		dev.off()
		pdf(paste0("Figure_H_NEW2.pdf"), width=11.0, height=4.5) # dev.new(width=11.0, height=4.5)
		par(mfrow=c(2,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		minVals1 = min(diag(matrices_mean2[[1]])); maxVals1 = max(diag(matrices_mean2[[1]]))
		mat = matrices_mean2[[1]]; diag(mat) = NA; minVals2 = min(mat, na.rm=T); maxVals2 = max(mat, na.rm=T)
		vec = branches_mean2[[1]]; minVals3 = min(vec, na.rm=T); maxVals3 = max(vec, na.rm=T)
		for (n in 2:nberOfReplicates)
			{
				mat1 = matrices_mean2[[n]]; mat2 = mat1
				diag(mat2) = NA; mat3 = branches_mean2[[n]]
				if (minVals1 > min(diag(mat1))) minVals1 = min(diag(mat1))
				if (maxVals1 < max(diag(mat1))) maxVals1 = max(diag(mat1))
				if (minVals2 > min(mat2,na.rm=T)) minVals2 = min(mat2,na.rm=T)
				if (maxVals2 < max(mat2,na.rm=T)) maxVals2 = max(mat2,na.rm=T)
				if (minVals3 > min(mat3,na.rm=T)) minVals3 = min(mat3,na.rm=T)
				if (maxVals3 < max(mat3,na.rm=T)) maxVals3 = max(mat3,na.rm=T)
			}
		for (n in 1:nberOfReplicates)
			{
				mat = matrices_mean2[[n]]; vec = branches_mean2[[n]]
				plot(NYboroughs, col=col_boroughs, border="gray60", lwd=0.5)
				points(centroids, cex=10*((diag(mat)-minVals1)/(maxVals1-minVals1)), pch=16, col="#4D4D4D50")
				if (plottingBoroughCladeBranches == TRUE)
					{
						points(centroids, cex=10*((vec-minVals1)/(maxVals1-minVals1)), pch=1, col="#4D4D4D75", lwd=0.5, lty=2)
					}
				for (i in 1:dim(NYboroughs)[1])
					{
						for (j in 1:dim(NYboroughs)[1])
							{
								if ((i != j)&(mat[i,j]>=1))
									{
										LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(mat[i,j]/maxVals2))+0.04
										curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
													lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
									}
							}
					}
				if (n == 1)
					{
						vS = 5; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.900), cbind(-74.12,40.900), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 10; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.875), cbind(-74.12,40.875), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						vS = 20; LWD = (((vS-minVals2)/(maxVals2-minVals2))*2)+0.1; arrow = (0.1*(vS/maxVals2))+0.04
						curvedarrow(cbind(-74.20,40.850), cbind(-74.12,40.850), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
									lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
						mtext("5", at=-74.11, line=-1.90, adj=0, col="gray30", cex=0.7)
						mtext("10", at=-74.11, line=-2.80, adj=0, col="gray30", cex=0.7)
						mtext("20", at=-74.11, line=-3.70, adj=0, col="gray30", cex=0.7)
						points(cbind(rep(-74.15,4),rep(40.75,4)), cex=10*((seq(5,20,5)-minVals1)/(maxVals1-minVals1)), pch=1, col="#4D4D4D", lwd=0.3)
					}
			}
		dev.off()
	}

	# B9. Visualising the spatio-temporal distribution of spike mutations

mutations = c("L5F","S98F","Y145H","D215Y/H","H1101Y") # 5 spike mutations that occur at least 3 times in the cohort (without D614G)
mutations_tab = read.csv("NY_mutations_file_2.csv"); metadata = read.csv("A_on_IQTREE-TreeTime/Sampling_NY_state.csv", head=T)
if (showingPlots)
	{
		wd = getwd(); nberOfReplicates = 10; mccs = list()
		setwd(paste0(wd,"/B_integrated_analyses/Replicate_RRW_analyses/"))
		for (n in 1:nberOfReplicates) mccs[[n]] = read.csv(paste0("Replicate_",i,"_MCC.csv"), head=T)
		minYear = min(mccs[[1]][,"startYear"]); maxYear = max(mccs[[1]][,"endYear"])
		for (n in 2:nberOfReplicates)
			{
				if (minYear > min(mccs[[n]][,"startYear"])) minYear = min(mccs[[n]][,"startYear"])
				if (maxYear < max(mccs[[n]][,"endYear"])) maxYear = max(mccs[[n]][,"endYear"])
			}
		colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(141)[16:116]); setwd(wd)
		pdf(paste0("Figure_J_NEW.pdf"), width=11.0, height=4.5) # dev.new(width=11.0, height=4.5)
		par(mfrow=c(2,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); cexNode = 1.2
		for (n in 1:length(mutations))
			{
				if (!grepl("\\/",mutations[n]))
					{
						seqIDs = mutations_tab[which(grepl(mutations[n],mutations_tab[,"spikeMutations"])),"seqName"]
					}	else	{
						seqIDs = mutations_tab[which(grepl("D215Y",mutations_tab[,"spikeMutations"])),"seqName"]
						seqIDs = c(seqIDs, mutations_tab[which(grepl("D215H",mutations_tab[,"spikeMutations"])),"seqName"])
					}
				sub = metadata[which(metadata[,"sequenceID"]%in%seqIDs),]
				endYears_indices = (((sub[,"collectionDate"]-minYear)/(maxYear-minYear))*100)+1
				endYears_colours = colourScale[endYears_indices]
				plot(NYboroughs, col="gray90", border="white", lwd=0.5)
				for (i in 1:dim(sub)[1])
					{
						points(sub[i,"longitude"], sub[i,"latitude"], pch=16, col=endYears_colours[i], cex=cexNode)
						points(sub[i,"longitude"], sub[i,"latitude"], pch=1, col="gray30", cex=cexNode, lwd=0.2)
					}
				mtext(mutations[n], side=1, at=-74.12, line=-10, col="gray30", cex=0.7)
			}
		dev.off()
	}

