# ##' @param f `character(1)` with the path to an MassBank file.
# ##'
# ##' @param msLevel `numeric(1)` with the MS level. Default is 2.
# ##'
# ##' @param metaDataBlocks `data.frame` data frame indicating which metadata shall
# ##'     be read
# ##'
# ##' @param nonStop `logical(1)` whether import should be stopped if an
# ##'     Massbank file does not contain all required fields. Defaults to
# ##'     `nonStop = FALSE`.
# ##'
# ##' @param ... Additional parameters, currently ignored.
# ##'
# ##' @importFrom S4Vectors DataFrame cbind.DataFrame
# ##' @importFrom IRanges NumericList
# ##'
# ##' @author Michael Witting
# ##'
# ##' @noRd
# .read_sirius <- function(f, msLevel = 2L,
#                            metaBlocks = metaDataBlocks(),
#                            nonStop = FALSE,
#                            ...) {
#
#   if (length(f) != 1L)
#     stop("Please provide a single Massbank file.")
#
#   mb <- scan(file = f, what = "",
#              sep = "\n", quote = "",
#              allowEscapes = FALSE,
#              quiet = TRUE)
#
#   begin <- grep(">compound", mb) + 1L
#   end <- grep("^//$", mb)
#
#   n <- length(begin)
#   spec <- vector("list", length = n)
#
#   for (i in seq(along = spec)) {
#     spec[[i]] <- .extract_mb_spectrum(mb[begin[i]:end[i]])
#   }
#
#   res <- DataFrame(do.call(rbind, spec))
#
#   for (i in seq_along(res)) {
#     if (all(lengths(res[[i]]) == 1))
#       res[[i]] <- unlist(res[[i]])
#   }
#
#   res$mz <- IRanges::NumericList(res$mz)
#   res$intensity <- IRanges::NumericList(res$intensity)
#   res$dataOrigin <- f
#   res$msLevel <- as.integer(msLevel)
#   res
#
# }
#
# ##' @param mb `character()` of lines defining a spectrum in mgf
# ##'     format.
# ##'
# ##' @importFrom utils tail type.convert
# ##'
# ##' @author Michael Witting
# ##'
# ##' @noRd
# .extract_sirius_ms1_spectrum <- function(mb, nonStop = FALSE) {
#
#   #read the spectrum
#   spectrum_start <- grep("PK$PEAK:", mb, fixed = TRUE) + 1
#   spectrum_end <- tail(grep("//", mb, fixed = TRUE), 1) - 1
#
#   splitted <- strsplit(mb[spectrum_start:(spectrum_end)]," ")
#   spectrum <- matrix(nrow = spectrum_end + 1 - spectrum_start, ncol = 3)
#
#   for(k in 1:length(splitted)){
#     splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
#     spectrum[k,] <- splitted[[k]]
#   }
#
#   # convert to data frame and adjust data type
#   spectrum <- as.data.frame(spectrum, stringsAsFactors = FALSE)
#   spectrum[] <- lapply(spectrum, type.convert)
#   colnames(spectrum) <- c("mz", "intensity", "rel.intensity")
#
#   # isolate Accession, name etc...
#   meta <- list()
#
#   meta$accession <- substring(grep("ACCESSION:", mb, value = TRUE, fixed = TRUE), 12)
#   meta$name <- as.list(substring(grep("CH$NAME:", mb, value = TRUE, fixed = TRUE), 10))
#   meta$smiles <- substring(grep("CH$SMILES:", mb, value = TRUE, fixed = TRUE), 12)
#   meta$exact_mass <- as.numeric(substring(grep("CH$EXACT_MASS:", mb, value = TRUE, fixed = TRUE), 16))
#   meta$formula <- substring(grep("CH$FORMULA:", mb, value = TRUE, fixed = TRUE), 13)
#   meta$iupac <- substring(grep("CH$IUPAC:", mb, value = TRUE, fixed = TRUE), 11)
#   meta$link_cas <- substring(grep("CH$LINK: CAS", mb, value = TRUE, fixed = TRUE), 14)
#   meta$link_pubchem <- substring(grep("CH$LINK: PUBCHEM", mb, value = TRUE, fixed = TRUE), 18)
#   meta$link_inchikey <- substring(grep("CH$LINK: INCHIKEY", mb, value = TRUE, fixed = TRUE), 19)
#   meta$link_chemspider <- substring(grep("CH$LINK: CHEMSPIDER", mb, value = TRUE, fixed = TRUE), 21)
#   meta$ms_col_energy <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_ENERGY", mb, value = TRUE, fixed = TRUE), 40)
#   meta$ms_frag_mode <- substring(grep("AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE", mb, value = TRUE, fixed = TRUE), 42)
#   meta$chrom_column <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_NAME", mb, value = TRUE, fixed = TRUE), 32)
#   meta$focus_precursor_type <- substring(grep("MS$FOCUSED_ION: PRECURSOR_TYPE", mb, value = TRUE, fixed = TRUE), 32)
#   meta$rtime_string <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME", mb, value = TRUE, fixed = TRUE), 35)
#
#   meta <- .cleanParsing(meta)
#
#   # type conversion
#   meta$ms_col_energy <- as.numeric(meta$ms_col_energy)
#
#   if(!is.na(meta$rtime_string)) {
#
#     rtime <- as.numeric(regmatches(meta$rtime_string, regexpr("[[:digit:]]+\\.[[:digit:]]+", meta$rtime_string)))
#     if(grepl("min", meta$rtime_string)) rtime <- rtime * 60
#
#   } else {
#
#     rtime <- NA_real_
#
#   }
#
#   precursorMz <- as.numeric(substring(grep("MS$FOCUSED_ION: PRECURSOR_M/Z", mb, value = TRUE, fixed = TRUE), 30))
#   precursorIntensity <- as.numeric(substring(grep("MS$FOCUSED_ION: PRECURSOR_INT", mb, value = TRUE, fixed = TRUE), 31))
#
#   # back up if no values are supplied
#   #if(!length(rtime)) rtime <- NA_real_
#   if(!length(precursorMz)) precursorMz <- NA_real_
#   if(!length(precursorIntensity)) precursorIntensity <- 100
#
#   title <- substring(grep("RECORD_TITLE:",
#                           mb,
#                           value = TRUE,
#                           fixed = TRUE),
#                      15)
#
#   list(accession = meta$accession,
#        name = meta$name,
#        smiles = meta$smiles,
#        exact_mass = meta$exact_mass,
#        formula = meta$formula,
#        iupac = meta$iupac,
#        link_cas = meta$link_cas,
#        link_pubchem = meta$link_pubchem,
#        link_inchikey = meta$link_inchikey,
#        link_chemspider = meta$link_chemspider,
#        ms_frag_mode = meta$ms_frag_mode,
#        chrom_column = meta$chrom_column,
#        focus_precursor_type = meta$focus_precursor_type,
#        rtime = rtime,
#        scanIndex = as.integer(1),
#        precursorMz = precursorMz,
#        precursorIntensity = precursorIntensity,
#        precursorCharge = as.integer(0),
#        mz = spectrum$mz,
#        intensity = spectrum$intensity,
#        collisionEnergy = meta$ms_col_energy,
#        title = title)
#
# }
#
#
# ##' Clean parsing
# ##'
# ##' @title Cleaning meta data parsing
# ##'
# ##' @param x `List` with entries
# ##'
# ##' @return `List` with cleaned entries
# ##'
# ##' @noRd
# .cleanParsing <- function(x) {
#
#   for(i in 1:length(x)) {
#
#     if(!length(x[[i]])) {
#       x[[i]] <- NA_character_
#     }
#
#   }
#
#   x
#
# }

##' write Sirius files
##'
##' @title Write Sirius files
##'
##'
export_sirius <- function(spectra,
                          f,
                          path) {

}

##' export single file
##'
##'
.export_sirius_single <- function(spectra,
                                  path) {

  # isolate MS1 and MS2 spectra
  ms1_spectra <- spectra[which(spectra$msLevel == 1L)]
  ms2_spectra <- spectra[which(spectra$msLevel == 2L)]

  # combine MS1 spectra in case there are multiple


}
