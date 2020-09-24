const bodyparser = require("body-parser");
const express = require("express");
const path = require("path");
const app = express();

const { spawn } = require("child_process");

const port = process.env.PORT || 3000;

app.set("views", path.join(__dirname));
app.set("view engine", "ejs");
app.use(bodyparser.json());
app.use(bodyparser.urlencoded({ extended: false }));
app.get("/", function (req, res) {
  res.render("sampleForm");
});
app.get("/about", function (req, res) {
  res.render("about");
});
app.post("/saveData", (req, res) => {
  var dataToSend;
  var final;
  var motif35 = JSON.parse("[" + req.body.seq35 + "]");
  var motif10 = JSON.parse("[" + req.body.seq10 + "]");

  const python = spawn("python", ["test.py", motif35, motif10, req.body.seq]);

  // collect data from script
  python.stdout.on("data", function (data) {
    console.log("Pipe data from python script ...");
    dataToSend = data.toString();
    final = eval(dataToSend);
  });
  // in close event we are sure that stream from child process is closed
  python.on("close", (code) => {
    console.log(`child process close all stdio with code ${code}`);
    res.status(200).render("results", { data: final });
  });
});

app.listen(port, () =>
  console.log(`Example app listening on port 
${port}!`)
);
