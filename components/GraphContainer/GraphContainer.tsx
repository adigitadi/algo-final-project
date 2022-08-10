import React, { useMemo, useState } from "react";
import CytoscapeComponent from "react-cytoscapejs";
import { Config, names, uniqueNamesGenerator } from "unique-names-generator";
import { v4 as uuidv4 } from "uuid";
import { Rider, Taxi, Variant } from "../../types";

import { createTaxiGraph, TaxiGraph } from "../../utils";

import RiderList from "../RiderList";

const config: Config = {
  dictionaries: [names],
};

const GraphContainer = () => {
  const [graphData, setGraphData] = useState<TaxiGraph>([]);
  const [riders, setRiders] = useState<Rider[]>([]);
  const [taxis, setTaxis] = useState<Taxi[]>([]);
  const [numTaxis, setNumTaxis] = useState(3);
  const [height, setHeight] = useState(12);
  const [width, setWidth] = useState(12);
  const [numRiders, setNumRiders] = useState(10);

  // HERE. Probably an array of string ids, in order (?)
  // if multiple paths, then that's going to be an array of arrays.
  const [path, setPath] = useState([]);

  // TODO: https://github.com/plotly/react-cytoscapejs/issues/46s

  // populates the inital graph, and any changes
  useMemo(() => {
    setRiders([]);
    setTaxis([]);
    setGraphData(createTaxiGraph(height, width));
  }, [height, width]);

  const calculatePaths = () => {
    console.log("NULL FUNCTION");
  };

  const populateRiders = () => {
    const riders: Rider[] = [];
    for (let i = 0; i < numRiders; i++) {
      const newRider = createRandomRider();
      riders.push(newRider);
    }
    setRiders(riders);
  };

  const createRandomRider = () => {
    const numNodes = height * width;
    const source = Math.floor(Math.random() * numNodes).toString();
    const characterName: string = uniqueNamesGenerator(config);
    const color = "#b0c4de";
    const destination = Math.floor(Math.random() * numNodes).toString();
    return {
      id: uuidv4(),
      name: characterName,
      source: source,
      destination: destination,
      color: color,
      variant: Variant.rider,
    };
  };

  const createRandomTaxi = () => {
    const source = 0;
    const characterName: string = uniqueNamesGenerator(config);
    const color = "#ffff00";
    return {
      id: uuidv4(),
      name: characterName,
      source: source.toString(),
      color: color,
      variant: Variant.taxi,
      capacity: 4,
      carrying: 0,
    };
  };

  const populateTaxis = () => {
    const taxis: Taxi[] = [];
    for (let i = 0; i < numTaxis; i++) {
      const newTaxi = createRandomTaxi();
      taxis.push(newTaxi);
    }
    setTaxis(taxis);
  };

  const clearRiders = () => {
    setRiders([]);
  };

  const clearTaxis = () => {
    setTaxis([]);
  };

  const onHeightChange = (event: { target: { value: any } }) => {
    setHeight(event.target.value);
  };

  const onWidthChange = (event: { target: { value: any } }) => {
    setWidth(event.target.value);
  };

  const onNumRidersChange = (event: { target: { value: any } }) => {
    setNumRiders(event.target.value);
  };

  const onNumTaxisChange = (event: { target: { value: any } }) => {
    setNumTaxis(event.target.value);
  };

  const cyRef = React.useRef<cytoscape.Core | undefined>();

  if (typeof window !== "undefined") {
    return (
      <div className="static  h-screen max-h-screen grid grid-cols-5">
        <RiderList
          height={height}
          width={width}
          riders={riders}
          numRiders={numRiders}
          taxis={taxis}
          numTaxis={numTaxis}
          onHeightChange={onHeightChange}
          onWidthChange={onWidthChange}
          onNumRidersChange={onNumRidersChange}
          onNumTaxisChange={onNumTaxisChange}
          populateRiders={populateRiders}
          clearRiders={clearRiders}
          populateTaxis={populateTaxis}
          clearTaxis={clearTaxis}
          calculatePaths={calculatePaths}
        />
        <CytoscapeComponent
          cy={(cy) => (cyRef.current = cy)}
          className="col-span-4"
          elements={graphData}
          stylesheet={[
            {
              selector: "node",
              style: {
                width: 10,
                height: 10,
                shape: "rectangle",
              },
            },
            {
              selector: "edge",
              style: {
                width: 1.5,
                "target-arrow-color": "red",
                "target-arrow-shape": "triangle",
                "curve-style": "taxi",
              },
            },
          ]}
          style={{ margin: "5px" }}
        />
      </div>
    );
  } else {
    return <></>;
  }
};

export default GraphContainer;
