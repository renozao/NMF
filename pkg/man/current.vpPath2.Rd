\name{current.vpPath2}
\alias{current.vpPath2}
\title{Custom Access to Viewport Paths}
\usage{
  current.vpPath2()
}
\description{
  This function was defined to substitute current.vpPath
  because this latter calls grid.Call("L_currentViewport")
  which calls 'L_gridDirty' before effectively calling
  'L_currentViewport'. This starts a new page (not sure
  why), which is not ok when writing to a device that
  records each page (e.g. PDF or SVG)
}

