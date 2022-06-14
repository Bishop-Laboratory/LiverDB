# UI function to make ? button
helpButton <- function(message) {
  return(
    add_prompt(
      ui_element = span(HTML('<i class="fa fa-question-circle"></i>')),
      message = message, position = "right"
    )
  )
}

# Make headers
makeHeaders <- function(title, message, fs=1.3) {
  tagList(
    span(span(title, style=paste0("font-size: ", fs, "em;")), helpButton(message))
  )
}
