library(shiny)
library(data.table)
library(dplyr)
library(DT)
library(openxlsx)
library(readr)
library(reshape2)
library(rsconnect)
library(shinyalert)
library(shinymanager)
library(shinythemes)
library(stringr)
library(tidyverse)
library(writexl)

options(stringsAsFactors = FALSE)
options(shiny.maxRequestSize = 300*1024^2)
options(scipen = 999) 

#Set controls for ddPCR Setup
numCons <- 3
blankVal <- rep(NA, numCons)

ddpcrCons <- data.frame(
  Well = "",
  experiment = blankVal,
  ddPCR.thermocycling.program = rep("Optimized", numCons),
  Sample.num = blankVal,
  Sample.name = c("ddPCR Thresholding Control", "ddPCR Positive Control", "Negative Control"),
  Well.type = rep("Individual", numCons),
  dilution = c(1, 100000, NA),
  volume_uL = c(1, 2, NA),
  target.organism = rep("MTB", numCons),
  RT.batch = blankVal,
  type = rep("ETS1_MGB, 23S", numCons),
  PCR_probe = rep("New Zen.IowaBlack probe", numCons),
  comments = blankVal,
  ddPCR.performer = blankVal,
  Replicates = blankVal,
  check.names = FALSE
)

################################################################################
################################################################################
################################################################################

ui <- navbarPage(
  theme = shinytheme("cerulean"),
  title = "ddPCR Setup",
  tabPanel(
    useShinyalert(),
    title = "Freezerworks Template",
    sidebarLayout(
      sidebarPanel(
        fileInput("samp_lst", "Select the ddPCR Sample List",
                  multiple = FALSE,
                  accept = ".csv"),
        div(style = "margin-top: -17px"),
        fileInput("frw_exprt", "Select the Most Recent FW Export",
                  multiple = FALSE,
                  accept = ".xlsx"),
        div(style = "margin-top: -12px"),
        splitLayout(cellWidths = c("50%", "50%"), 
                    textInput("exp_name", label = NULL, placeholder = "Enter experiment name"), 
                    textInput("user_name", label = NULL, placeholder = "Enter ddPCR performer")
        ),
        downloadButton('downloadSetup', 'ddPCR Setup Files'),
        h6(""),
        radioButtons("frw_file_btn", 
                     "Choose File to Display",
                     choiceNames = c("ddPCR Inventory Merge",
                                     "Freezerworks Update",
                                     "ddPCR Master Index Template"),
                     choiceValues = c("fzrw_mrg",
                                      "fzrw_up",
                                      "ddpcr_mi")
        ),
        h6(""),
        strong("1. Inventory merge:"), 
        p("Combines sample list with the Freezerworks inventory"),
        strong("2. FW Update: "), 
        p("Sample list import file for Freezerworks updates"),
        strong("3. ddPCR Master Index Template:"), 
        p("Formatted sample list for use in ddPCR master index"),
        div(style = "margin-top: -20px"),
        width = 5
      ),
      mainPanel(DT::dataTableOutput("fzrwTable"), width = 7)
    )
  )
)

################################################################################
################################################################################
################################################################################

server <- function(input, output, session) {
  
  #Import ddPCR sample list 
  fzrw_inven <- reactive({
    req(input$samp_lst)
    req(input$frw_exprt)
    withProgress(message = 'Generating Files', value = 0, {
      incProgress(1)
      #Load the ddPCR sample list .csv file
      ddpcr_lst <- read.csv(input$samp_lst$datapath,
                            na.strings = c("", "NA"),
                            fileEncoding="latin1",
                            row.names = NULL)
      #Load the FW inventory export file
      fzrw_export <- read.xlsx(input$frw_exprt$datapath)
      
      ddpcr_lst <- ddpcr_lst %>%
        filter(!is.na(Sample.num)) %>%
        #Rename some of the columns for readability
        rename(
          Dilution = Adjusted.Recommend..Dil,
          Recommended.uL = Adjusted.Recommend..uL) %>%
        #Convert the dilutions to match the ones used in FW
        mutate(
          #Create a new column that stores the dilution factor as a numerical
          #value. This will be used to sort the sample list later, as this cannot
          #be done when the dilutions are in string form as they will be in 
          #the "Dilution" column.
          Dilution = gsub(",", "", Dilution),
          DilutionNum = as.numeric(Dilution),
          Dilution = as.character(paste0("1:", Dilution)),
          Dilution = ifelse(Dilution == "1:10000000000", "1:10B", Dilution),
          Dilution = ifelse(Dilution == "1:1000000000", "1:1B", Dilution),
          Dilution = ifelse(Dilution == "1:100000000", "1:100M", Dilution),
          Dilution = ifelse(Dilution == "1:10000000", "1:10M", Dilution),
          Dilution = ifelse(Dilution == "1:1000000", "1:1M", Dilution),
          Dilution = ifelse(Dilution == "1:100000", "1:100K", Dilution),
        )
      #If a replicates column is not included in sample list, replicates assumed to = 1
      if(!"Replicates" %in% colnames(ddpcr_lst)){
        ddpcr_lst$Replicates <- 1
      }

      incProgress(1)
      
      #Format the inventory export
      fzrw_export <- fzrw_export %>% 
        select(
          Sample.num,
          Aliquot.type,
          `RT.-.RT.Batch`,
          Dilution,
          Number.of.Thaws,
          Current.Amount,
          Freezer.Section,
          Position.1,
          Position.2,
          Position.3,
          Position.4,
          Position.5,
          Unique.Aliquot.ID
        ) %>%
        filter(Aliquot.type == "cDNA" & !is.na(Dilution)) %>%
        rename(
          RT.batch = `RT.-.RT.Batch`
        ) %>%
        mutate(
          #Change dilution factor == "neat" to "1:1" so that it matches ddPCR sample list
          Dilution = ifelse(Dilution == "Neat", "1:1", Dilution),
          RT.batch = as.integer(RT.batch),
          Sample.num = as.integer(Sample.num)
        ) %>%
        select(-Aliquot.type)
      
      #merge the ddPCR sample list and FW inventory export
      fzrw_merged <- merge(ddpcr_lst, fzrw_export, by = c("Sample.num", "Dilution", "RT.batch"), all.x = TRUE)
      
      incProgress(1)
      
      #Sort the sample list/FW merge by ddPCR recommended uL, dilution, then Sample.num,
      #this is how the assays are organized for the experiment. 
      fzrw_merged <- fzrw_merged %>%
        arrange(
          Recommended.uL,
          DilutionNum,
          Sample.num
        ) %>%
        relocate(
          RT.batch, .after = Sample.name
        ) %>%
        relocate(
          Dilution, .after = RT.batch
        )
      
      #Generate a list of samples that failed to merge for some reason. There
      #are three possible reasons this may occur that are accounted for:
      # 1. The RT batch provided in the sample list cannot be associated with the
      # Sample.num in the FW inventory export
      # 2. The necessary dilution associated with the Sample.num/RT batch combination 
      # in the assay does not exist in the FW inventory export.
      # 3. There is currently not enough cDNA to perform the assay at the 
      # recommended dilution/volume. Technically, this doesn't cause the merge
      # to fail, but we'll treat it as if it does to alert the user to the problem
      
      # When the merge fails, the freezer position is NA, so use that to find
      # "missing" samples. Then, we'll add to that list any samples that have
      # insufficient volume of cDNA for the assay
      missing_samples <- fzrw_merged %>%
        filter(
          is.na(Freezer.Section) |
            Current.Amount < ((Recommended.uL + (Recommended.uL * 0.25)) * Replicates)
        ) %>%
        select(Sample.num,
               RT.batch,
               Dilution) %>%
        arrange(Sample.num, RT.batch, Dilution)
      
      #Get all samples *from the FW inventory export* that match Sample.nums in
      #the missing_samples data. This is necessary to determine why the merge
      #failed. 
      found_samples <- fzrw_export %>%
        filter(Sample.num %in% missing_samples$Sample.num) %>%
        arrange(Sample.num, RT.batch, Dilution)
      
      incProgress(1)
      
      #initiate a list that will hold Sample.nums in the missing_samples data, as well
      #as a description of why those samples were not found (i.e. the necessary
      #dilution does not exist)
      failure_msgs = list()
      
      # Now, we're going to iterate through all the missing samples, determine
      # why they are in the missing sample list, and then add those Sample.nums
      # and the reason for the merge failure to the list above. If there are
      # no missing samples, this step is skipped. 
      if(!all(is.na(missing_samples$Sample.num))){
        for (row in 1:nrow(missing_samples)) {
          #First, we'll check if the RT batch is missing
          checkBatch = missing_samples[row, "RT.batch"]
          temp_data <- found_samples %>%
            filter(Sample.num == missing_samples[row, "Sample.num"] & 
                     RT.batch == missing_samples[row, "RT.batch"])
          if(!checkBatch %in% temp_data$RT.batch){
            msg <- paste0( 
              missing_samples[row, "Sample.num"],
              " : RT batch\n")
            failure_msgs <- append(failure_msgs, msg)
            missing_samples[row, "Failure"] <- "RT.Batch Not Found"
            #If the RT batch is found, then we'll check if the dilution is missing
          } else if(!missing_samples[row, "Dilution"] %in% temp_data$Dilution) {
            msg <- paste0( 
              missing_samples[row, "Sample.num"],
              " : Dilution\n")
            failure_msgs <- append(failure_msgs, msg)
            missing_samples[row, "Failure"] <- "Dilution Not Found"
            #Finally, we'll check if there is sufficient volume to run the experiment
          } else {
            msg <- paste0( 
              missing_samples[row, "Sample.num"],
              ": Volume\n")
            failure_msgs <- append(failure_msgs, msg)
            missing_samples[row, "Failure"] <- "Insufficient Volume"
          }
        }
      }
      # If there are ay samples in the "missing_samples" dataframe, we'll add
      # The "failure" column that is produced above to the merged dataframe
      # so that the user can see why the merge failed, and we'll also create
      # a Shiny alert to make this more obvious. 
      if(!all(is.na(missing_samples$Sample.num))){
        fzrw_merged <- merge(fzrw_merged, 
                             missing_samples[, c("Sample.num", "Failure")], by = "Sample.num", all.x = TRUE)
        
        shinyalert(paste("Merge failed for samples:", paste(failure_msgs, collapse="")), type = "warning")
        
      } else{
        fzrw_merged$Failure <- NA
      }
      #Rename the columns from the inventory dataframe to make the more user-friendly
      fzrw_merged <- fzrw_merged %>%
        mutate(Failure = ifelse(is.na(Failure), "Found", Failure)) %>%
        rename(
          Sample.Status = Failure,
          Freezer = Freezer.Section,
          Shelf = Position.1,
          Rack = Position.2,
          Box = Position.3,
          Row = Position.4,
          Column = Position.5
        ) %>%
        #arrange the samples based on how they are loaded onto the ddPCR plate
        arrange(
          Recommended.uL,
          DilutionNum,
          Sample.num
        ) %>%
        #filter to relevant columns
        select(
          Sample.Status,
          Sample.num,
          Sample.name,
          RT.batch,
          Dilution,
          Recommended.uL,
          Replicates,
          Current.Amount,
          Freezer,
          Shelf,
          Rack,
          Box,
          Row,
          Column,
          Unique.Aliquot.ID,
          Number.of.Thaws
        ) 
    })
  })
  
  # Now, we have a dataframe that contains the original sample list merged with
  # the inventory export file so that the user can easily find the samples. One
  # remaining challenge is that the current volume of the sample and the number
  # of times that sample needs to be thawed must be accounted for in the 
  # inventory software. In order to avoid needing to find all the samples and
  # update them individually (a painstaking and error-prone process), we'll
  # leverage the samples unique ID, which is automatically associated with
  # a sample by the inventory software, to create a dataframe of the samples'
  # unique ID as well as the current volume and number of thaws, which are
  # automatically updated based on the conditions of the experiment.
  fzrw_update <- reactive({
    req(fzrw_inven())
    fzrw_upd <- fzrw_inven() %>%
      mutate(
        Current.Amount = Current.Amount - ((Recommended.uL + (Recommended.uL * 0.25)) * Replicates),
        Current.Amount = ifelse(Current.Amount < 0 , "Insufficient Volume", Current.Amount),
        Unique.Aliquot.ID = ifelse(is.na(Unique.Aliquot.ID), 
                                   paste0("Missing Sample.num: ", Sample.num), 
                                   Unique.Aliquot.ID),
        Number.of.Thaws = Number.of.Thaws + 1
      ) %>%
      select(
        Unique.Aliquot.ID,
        Current.Amount,
        Number.of.Thaws
      ) %>%
      rename(
        `Number of Thaws` = Number.of.Thaws
      )
  })
  
  # The last step in the ddPCR setup is to create a file that properly formats
  # the sample list and associated metadata for the master index that contains
  # information of all samples that have ever been run using this procedure. 
  # This is a crucial step for integrating experiment data and sample metadata. 
  # Previously, this was done manually by copy/pasting and text entry, which is 
  # a cumbersome and error-prone process. 
  
  ddpcr_MI <- reactive({
    req(fzrw_inven())
    ddpcr <- fzrw_inven() %>%
      #Create and populate columns to match those that exist in the master index
      mutate(
        Well = "",
        experiment = as.character(input$exp_name),
        ddPCR.thermocycling.program = "Optimized",
        Well.type = "Individual",
        target.organism = "MTB",
        type = "ETS1_MGB, 23S", 
        PCR_probe = "New Zen.IowaBlack probe",
        comments = "",
        ddPCR.performer = input$user_name
      ) %>%
      rename(
        dilution = Dilution,
        volume_uL = Recommended.uL
      ) %>%
      select(
        Well,
        experiment,
        ddPCR.thermocycling.program,
        Sample.num,
        Sample.name,
        Well.type,
        dilution,
        volume_uL,
        target.organism,
        RT.batch,
        type,
        PCR_probe,
        comments,
        ddPCR.performer,
        Replicates
      )
    
    #Initialize an empty dataframe that contains the same column names as above
    
    mastIdxTemp<- ddpcr[0,]
    
    # One challenge here is that samples are sometimes run in replicates. That 
    # needs to be reflected in the master index as duplicate rows that are repeated
    # as many times as there are replicates for that sample in the experiment. Once
    # the replicate rows are added, we need to add one last row indicating that
    # data from the above rows needs to be "merged".
    
    for(i in 1:nrow(ddpcr)) {
      numReps <- ddpcr[i, "Replicates"]
      row <- ddpcr[i, ]
      for(j in 1:numReps){
        mastIdxTemp <- rbind(mastIdxTemp, row)
      }
      if (numReps > 1){
        row$Well.type <- "Merge"
        mastIdxTemp <- rbind(mastIdxTemp, row)
      }
    }
    
    #Add rows for the control values 
    mastIdxTemp <- rbind(mastIdxTemp, ddpcrCons)
    
    # Filter to the row containing information for the negative controls. 
    # We will repeat this a certain number of times to "fill" the plate 
    # The plate is not technically full, but all wells in a column of the
    # plate must be full, which means that the total of the rows (***excluding
    # rows indicating the "merge"***) must be divisible by 8. We also need
    # to account for the fact that at least two negative controls must be 
    # included in each experiment. 
    Negs <- ddpcrCons %>%
      filter(Sample.name == "Negative Control")
    
    while(sum(str_count(mastIdxTemp$Well.type, "Individual")) %% 8 != 0) {
      mastIdxTemp <- rbind(mastIdxTemp, Negs)
    }
    
    if(sum(str_count(mastIdxTemp$Sample.name, "Negative Control")) < 2) {
      mastIdxTemp <- rbind(mastIdxTemp, Negs)
      while(sum(str_count(mastIdxTemp$Well.type, "Individual")) %% 8 != 0) {
        mastIdxTemp <- rbind(mastIdxTemp, Negs)
      }
    }
    
    wellIdx = 1
    lettIdx = 1
    
    # Now, we'll populate the column that indicates the plate well for a given 
    # sample, which has the format "A01", "B01", "C01"... These are only included
    # for non-merge indicating rows
    for (idx in 1:nrow(mastIdxTemp)){
      if(mastIdxTemp[idx, "Well.type"] == "Individual"){
        if(lettIdx <= 8){
          wellNum <- paste0("0", as.character(wellIdx))
          wellVal <- paste0(LETTERS[lettIdx], wellNum)
          mastIdxTemp[idx, "Well"] <- wellVal
          lettIdx <- lettIdx + 1
        } else{
          lettIdx <- 1
          wellIdx <- wellIdx + 1
          wellNum <- paste0("0", as.character(wellIdx))
          wellVal <- paste0(LETTERS[lettIdx], wellNum)
          mastIdxTemp[idx, "Well"] <- wellVal
          lettIdx <- lettIdx + 1
        }
      }
    }
    
    #Lastly, we'll convert the dilution factor to a numerical value so that
    #math can be performed on the samples based on that dilution factor. 
    mastIdxTemp <- mastIdxTemp %>% 
      select(-Replicates) %>%
      mutate(
        dilution = ifelse(dilution == "1:1", 1, dilution),
        dilution = ifelse(dilution == "1:10", 10, dilution),
        dilution = ifelse(dilution == "1:100", 100, dilution),
        dilution = ifelse(dilution == "1:1000", 1000, dilution),
        dilution = ifelse(dilution == "1:10000", 10000, dilution),
        dilution = ifelse(dilution == "1:100K", 100000, dilution),
        dilution = ifelse(dilution == "1:1M", 1000000, dilution),
        dilution = ifelse(dilution == "1:10M", 10000000, dilution),
        dilution = ifelse(dilution == "1:100M", 100000000, dilution),
        dilution = ifelse(dilution == "1:1B", 1000000000, dilution),
        dilution = ifelse(dilution == "1:10B", 10000000000, dilution)
      )
    return(mastIdxTemp)
  })
  
  ########################################################################
  ########################################################################
  ########################################################################
  #Display processed file
  
  output$firstData <- DT::renderDataTable(datatable(createTable(),
                                                    extensions = "FixedColumns", 
                                                    options = list(rownames= FALSE)
  ))
  
  output$fzrwTable <- DT::renderDataTable(
    
    if(input$frw_file_btn == "fzrw_up") {
      return(
        datatable(fzrw_update(),
                  options = list(
                    rownames= FALSE, 
                    paging = FALSE,
                    scrollY = "475px",
                    columnDefs = list(list(className = 'dt-center', targets = "_all"))
                  )
                  #rownames= FALSE,
        )
      )
    } else if (input$frw_file_btn == "fzrw_mrg") {
      return(
        datatable(fzrw_inven(),
                  options = list( 
                    scrollY= "475px",
                    scrollX = TRUE,
                    paging = FALSE,
                    autoWidth = TRUE,
                    searching = FALSE,
                    columnDefs = list(
                      list(className = 'dt-center', targets = "_all")
                    )
                  ), 
                  rownames = FALSE
        )  %>% formatStyle(
          'Replicates',
          target = "row",
          backgroundColor = styleEqual(c(2:6), 
                                       c(rep('#DBF3FA', 5)))) %>% 
          formatStyle(
            'Sample.Status',
            target = "row",
            backgroundColor = styleEqual(c("Dilution Not Found", "RT.Batch Not Found", "Insufficient Volume"), 
                                         c('#FFDFBF', '#FFDFBF', '#FFDFBF'))) %>%
          formatStyle("Sample.Status","white-space"="nowrap")
      )
    } else if(input$frw_file_btn == "ddpcr_mi") {
      return(
        datatable(ddpcr_MI()
                  # options = list(
                  #   paging = FALSE,
                  #   scrollY= "475px",
                  #   scrollX = TRUE)
                  # rownames= FALSE,
        )
      )
    }
  )
  
  ########################################################################
  ########################################################################
  ########################################################################
  #Download processed files
  output$downloadSetup <- downloadHandler(
    filename = function() {
      paste0(input$exp_name, " ddPCR Setup Files", ".zip")
    },
    content = function( file){
      
      # Set temporary working directory
      twd <- setwd(tempdir())
      on.exit(setwd(twd))
      
      # Save output files to the temporary directory
      writexl::write_xlsx(fzrw_inven(), paste0(input$exp_name, " - Inventory Merge.xlsx"))
      write.csv(fzrw_update(), paste0(input$exp_name, " - Freezerworks Update.csv"), row.names = FALSE)
      write.csv(ddpcr_MI(), paste0(input$exp_name, " - ddPCR Master Index Output.csv"), row.names = FALSE)
      
      # Stores outputs in Zip file
      zip(file, c(
        paste0(input$exp_name, " - Inventory Merge.xlsx"),
        paste0(input$exp_name, " - Freezerworks Update.csv"),
        paste0(input$exp_name, " - ddPCR Master Index Output.csv")
      ))
    }
  )
}

shinyApp(ui = ui, server = server)