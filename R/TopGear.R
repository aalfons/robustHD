# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Top Gear car data
#'
#' @description
#' The data set contains information on cars featured on the website of the
#' popular BBC television show \emph{Top Gear}.
#'
#' @usage
#' data("TopGear")
#'
#' @format
#' A data frame with 297 observations on the following 32 variables.
#' \describe{
#'   \item{\code{Maker}}{factor; the car maker.}
#'   \item{\code{Model}}{factor; the car model.}
#'   \item{\code{Type}}{factor; the exact model type.}
#'   \item{\code{Fuel}}{factor; the type of fuel (\code{"Diesel"} or
#'   \code{"Petrol"}).}
#'   \item{\code{Price}}{numeric; the list price (in UK pounds)}
#'   \item{\code{Cylinders}}{numeric; the number of cylinders in the engine.}
#'   \item{\code{Displacement}}{numeric; the displacement of the engine (in
#'   cc).}
#'   \item{\code{DriveWheel}}{factor; the type of drive wheel (\code{"4WD"},
#'   \code{"Front"} or \code{"Rear"}).}
#'   \item{\code{BHP}}{numeric; the power of the engine (in bhp).}
#'   \item{\code{Torque}}{numeric; the torque of the engine (in lb/ft).}
#'   \item{\code{Acceleration}}{numeric; the time it takes the car to get from
#'   0 to 62 mph (in seconds).}
#'   \item{\code{TopSpeed}}{numeric; the car's top speed (in mph).}
#'   \item{\code{MPG}}{numeric; the combined fuel consuption (urban + extra
#'   urban; in miles per gallon).}
#'   \item{\code{Weight}}{numeric; the car's curb weight (in kg).}
#'   \item{\code{Length}}{numeric; the car's length (in mm).}
#'   \item{\code{Width}}{numeric; the car's width (in mm).}
#'   \item{\code{Height}}{numeric; the car's height (in mm).}
#'   \item{\code{AdaptiveHeadlights}}{factor; whether the car has adaptive
#'   headlights (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{AdjustableSteering}}{factor; whether the car has adjustable
#'   steering (\code{"no"} or \code{"standard"}).}
#'   \item{\code{AlarmSystem}}{factor; whether the car has an alarm system
#'   (\code{"no/optional"} or \code{"standard"}).}
#'   \item{\code{Automatic}}{factor; whether the car has an automatic
#'   transmission (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{Bluetooth}}{factor; whether the car has bluetooth
#'   (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{ClimateControl}}{factor; whether the car has climate control
#'   (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{CruiseControl}}{factor; whether the car has cruise control
#'   (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{ElectricSeats}}{factor; whether the car has electric seats
#'   (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{Leather}}{factor; whether the car has a leather interior
#'   (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{ParkingSensors}}{factor; whether the car has parking sensors
#'   (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{PowerSteering}}{factor; whether the car has power steering
#'   (\code{"no"} or \code{"standard"}).}
#'   \item{\code{SatNav}}{factor; whether the car has a satellite navigation
#'   system (\code{"no"}, \code{"optional"} or \code{"standard"}).}
#'   \item{\code{ESP}}{factor; whether the car has ESP (\code{"no"},
#'   \code{"optional"} or \code{"standard"}).}
#'   \item{\code{Verdict}}{numeric; review score between 1 (lowest) and 10
#'   (highest).}
#'   \item{\code{Origin}}{factor; the origin of the car maker (\code{"Asia"},
#'   \code{"Europe"} or \code{"USA"}).}
#' }
#'
#' @source
#' The data were scraped from \code{http://www.topgear.com/uk/} on 2014-02-24.
#' Variable \code{Origin} was added based on the car maker information.
#'
#' @examples
#' data("TopGear")
#' summary(TopGear)

"TopGear"
