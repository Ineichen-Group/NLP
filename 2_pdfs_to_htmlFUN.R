#Conversion of pdfs in html

pdf_source_dir <- "All_PDFs_automated_RoB"
html_output_dir <- "All_PDFs_automated_RoB/html"

#run 'pip install PyMuPDF' in the terminal needed
fitz <- import("fitz")

pdf_files <- list.files(pdf_source_dir, pattern = "*.pdf", full.names = TRUE)

for (file_path in pdf_files) {
  file_name_html <- gsub(".pdf", ".html", basename(file_path))
  doc <- fitz$open(file_path)
  doc_text <- '<div id="page0" style="width:595.3pt;height:793.7pt">'
  for (p_nr in 0:(doc$page_count - 1)) {
    page <- doc$load_page(p_nr)
    text <- page$get_text("html")
    text <- gsub("<div.*?>|</div>", "", text, perl = TRUE)
    text <- gsub("<img.*?>", "", text, perl = TRUE) # trying to get rid of images
    doc_text <- paste(doc_text, text, sep = "")
  }
  doc_text <- paste(doc_text, "</div>", sep = "")
  write(doc_text, file.path(html_output_dir, file_name_html))
  doc$close()
}