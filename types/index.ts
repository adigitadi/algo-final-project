export enum Variant {
  taxi = "taxi",
  rider = "rider",
}

type Unit = {
  id: string;
  name: string;
  source: string;
  color: string;
  variant: Variant;
};

export type Rider = Unit & {
  destination: string;
};

export type Taxi = Unit & {
  capacity: number;
  carrying: number;
};
